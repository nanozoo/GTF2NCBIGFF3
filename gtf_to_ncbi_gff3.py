#!/usr/bin/env python3
"""
gtf_to_ncbi_gff3.py

Convert a StringTie-like GTF to an NCBI-style GFF3.
- Builds a gene -> mRNA -> exon hierarchy from 'gene', 'transcript', and 'exon' entries.
- Preserves useful description fields by folding them into Note=, Dbxref=, product= where possible.
- If CDS entries exist in the GTF, they are preserved and attached to the proper mRNA with phase.
- If no CDS are present, outputs gene/mRNA/exon only (appropriate for non-coding or draft annotations).

Usage:
  python gtf_to_ncbi_gff3.py input.gtf > output.gff3
  # or
  python gtf_to_ncbi_gff3.py input.gtf -o output.gff3
"""

import sys, re, argparse, collections
from urllib.parse import quote

AttrDict = dict

def parse_gtf_attributes(attr_field: str) -> AttrDict:
    """
    Parse a GTF attributes column into a dict of {key: [values...]}
    GTF attributes are key "value"; pairs, but sometimes stray tokens appear.
    We'll capture key "value" pairs and also keep a raw string for free-text mining.
    """
    d = collections.defaultdict(list)
    for m in re.finditer(r'(\S+)\s+"([^"]*)"', attr_field):
        d[m.group(1)].append(m.group(2))
    d["_raw"] = [attr_field]
    return d

def gff3_escape(v: str) -> str:
    """
    GFF3 requires escaping of special characters in attribute values using percent-encoding
    for ; = , and whitespace outside of unreserved set. Use urllib.parse.quote with safe set.
    """
    # Keep common safe ASCII. Encode semicolon, comma, equals, and spaces if present.
    return quote(v, safe=":+._-")  # allow colon, plus, dot, underscore, dash

def join_attrs(pairs):
    """
    Build a GFF3 attribute string from (key, value) with proper escaping.
    Supports list values by joining with comma (escaped as needed).
    """
    parts = []
    for k, v in pairs:
        if v is None or v == "":
            continue
        if isinstance(v, (list, tuple)):
            vv = ",".join(gff3_escape(str(x)) for x in v if x is not None and str(x) != "")
        else:
            vv = gff3_escape(str(v))
        parts.append(f"{k}={vv}")
    return ";".join(parts)

def first(d: AttrDict, *keys, default=None):
    for k in keys:
        if k in d and d[k]:
            return d[k][0]
    return default

def collect_dbxrefs(raw_attrs: str):
    """
    Mine Pfam-like accessions from free text. Return list of Dbxref items (e.g., Pfam:PF03031.14).
    Also pass through any 'Dbxref=...' like substrings if present (rare in GTF).
    """
    xrefs = set()
    # Pfam accessions
    for acc in re.findall(r"(Pfam:PF\d+(?:\.\d+)?)", raw_attrs):
        xrefs.add(acc)
    # TIGRFAM, CDD, InterPro (if present)
    for acc in re.findall(r"(TIGRFAMs?:TIGR\d+)", raw_attrs):
        xrefs.add(acc.replace("TIGRFAMs:", "TIGRFAM:"))
    for acc in re.findall(r"(CDD:cd\d+)", raw_attrs, flags=re.IGNORECASE):
        xrefs.add(acc.upper().replace("CDD:CD", "CDD:cd"))
    for acc in re.findall(r"(InterPro:IPR\d+)", raw_attrs):
        xrefs.add(acc)
    # Existing Dbxref='X:Y' literal fragments
    for acc in re.findall(r"Dbxref=['\"]([^'\"]+)['\"]", raw_attrs):
        xrefs.add(acc)
    return sorted(xrefs)

def summarize_notes(attrs: AttrDict):
    """Fold recognisable description fields into a compact Note string."""
    note_bits = []
    for key in sorted(attrs.keys()):
        if key.startswith("_"):  # skip raw
            continue
        low = key.lower()
        if any(tok in low for tok in ["desc", "hmm", "match", "note"]):
            # concatenate unique values
            vals = list(collections.OrderedDict.fromkeys(attrs[key]))
            s = f"{key}=" + ",".join(vals[:3])
            # limit explosion
            if len(vals) > 3:
                s += f",(+{len(vals)-3} more)"
            note_bits.append(s)
    note = "; ".join(note_bits)[:900]  # keep Note manageable
    return note

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gtf", help="Input GTF")
    ap.add_argument("-o", "--out", help="Output GFF3 file (default: stdout)")
    ap.add_argument("--source", default="GTF2NCBI", help="Source field to use in output")
    ap.add_argument("--locus-tag-prefix", dest="lt_prefix", default=None, help="If set, generate locus_tag like PREFIX_00001, PREFIX_00002 instead of using gene_id")
    ap.add_argument("--default-product", default="hypothetical protein", help="Fallback product for CDS if none found")
    args = ap.parse_args()

    out = sys.stdout if not args.out else open(args.out, "w", encoding="utf-8")
    print("##gff-version 3", file=out)

    # Accumulate records by gene->transcript->features
    genes = {}  # gene_id -> dict(meta, rows)
    transcripts = collections.defaultdict(dict)  # (gene_id, transcript_id) -> meta

    with open(args.gtf, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, source, feature, start, end, score, strand, frame, attr_field = parts
            attrs = parse_gtf_attributes(attr_field)

            gid = first(attrs, "gene_id", "geneid", "ID")
            tid = first(attrs, "transcript_id", "transcriptid", "Parent")
            rid = first(attrs, "exon_id", "exonid", "ID")

            # gene record
            if feature.lower() == "gene":
                # store gene bounds and attributes
                genes[gid] = {
                    "seqid": seqid, "start": int(start), "end": int(end),
                    "score": score if score != "." else ".",
                    "strand": strand, "frame": ".",
                    "attrs": attrs,
                }
            elif feature.lower() in ("transcript", "mrna"):
                key = (gid, tid)
                transcripts[key] = {
                    "seqid": seqid, "start": int(start), "end": int(end),
                    "score": score if score != "." else ".",
                    "strand": strand, "frame": ".",
                    "attrs": attrs,
                    "exons": [], "cds": []
                }
            elif feature.lower() == "exon":
                key = (gid, tid)
                exon_number = first(attrs, "exon_number", "exonNumber")
                transcripts[key].setdefault("exons", []).append({
                    "seqid": seqid, "start": int(start), "end": int(end),
                    "score": score if score != "." else ".",
                    "strand": strand, "frame": ".",
                    "attrs": attrs, "exon_number": exon_number
                })
            elif feature.lower() == "cds":
                key = (gid, tid)
                phase = frame if frame in ("0","1","2") else "."
                transcripts[key].setdefault("cds", []).append({
                    "seqid": seqid, "start": int(start), "end": int(end),
                    "score": score if score != "." else ".",
                    "strand": strand, "frame": phase,
                    "attrs": attrs
                })
            else:
                # ignore other feature types
                pass

    # Emit
    def emit(seqid, source, type_, start, end, score, strand, phase, attrs_pairs):
        print("\t".join([
            str(seqid), args.source, type_, str(start), str(end),
            score, strand, phase, join_attrs(attrs_pairs)
        ]), file=out)

    counter = 0
    for gid, g in genes.items():
        seqid = g["seqid"]; strand = g["strand"]
        gene_id = f"gene-{gid}"
        counter += 1
        if args.lt_prefix:
            locus_tag = f"{args.lt_prefix}_{counter:05d}"
        else:
            locus_tag = gid  # use gene_id as locus_tag
        gene_name = first(g["attrs"], "gene_name", "gene_symbol", default=gid)

        raw = g["attrs"].get("_raw", [""])[0]
        dbx = collect_dbxrefs(raw)
        note = summarize_notes(g["attrs"])

        gene_attrs = [
            ("ID", gene_id),
            ("Name", gene_name),
            ("locus_tag", locus_tag),
            ("gbkey", "Gene"),
        ]
        if dbx: gene_attrs.append(("Dbxref", dbx))
        if note: gene_attrs.append(("Note", note))

        emit(seqid, args.source, "gene", g["start"], g["end"], g["score"], strand, ".", gene_attrs)

        # find transcripts for this gene_id
        for (tgid, tid), t in sorted(transcripts.items()):
            if tgid != gid:
                continue
            rna_id = f"rna-{tid}"
            # bounds: if not provided by transcript row (rare), compute from exons
            if not t.get("start") or not t.get("end"):
                if t["exons"]:
                    starts = [e["start"] for e in t["exons"]]
                    ends = [e["end"] for e in t["exons"]]
                    t["start"], t["end"] = min(starts), max(ends)

            traw = t["attrs"].get("_raw", [""])[0]
            dbx_t = sorted(set(dbx) | set(collect_dbxrefs(traw)))
            note_t = summarize_notes(t["attrs"])

            mrna_attrs = [
                ("ID", rna_id),
                ("Parent", gene_id),
                ("Name", tid),
                ("gbkey", "mRNA"),
                ("transcript_id", tid),
            ]
            if dbx_t: mrna_attrs.append(("Dbxref", dbx_t))
            if note_t: mrna_attrs.append(("Note", note_t))

            emit(t["seqid"], args.source, "mRNA", t["start"], t["end"], t["score"], t["strand"], ".", mrna_attrs)

            # Exons
            # Sort by genomic order respecting strand
            exons_sorted = sorted(t["exons"], key=lambda e: (e["start"], e["end"]), reverse=(t["strand"]=="-"))
            for i, e in enumerate(exons_sorted, start=1):
                ex_num = e.get("exon_number") or str(i)
                exon_id = f"exon-{tid}.{ex_num}"
                eattrs = [
                    ("ID", exon_id),
                    ("Parent", rna_id),
                    ("gbkey", "mRNA"),
                    ("exon_number", ex_num),
                ]
                emit(e["seqid"], args.source, "exon", e["start"], e["end"], e["score"], e["strand"], ".", eattrs)

            # CDS (if present)
            if t["cds"]:
                cds_sorted = sorted(t["cds"], key=lambda c: (c["start"], c["end"]), reverse=(t["strand"]=="-"))
                # product inference from attributes if any
                prod = first(t["attrs"], "product", "protein_product", default=None)
                if not prod:
                    # Try heuristics from HMM descriptions
                    # take first HMMER_*_desc text or "hypothetical protein"
                    cand = None
                    for k in t["attrs"]:
                        if "desc" in k.lower() and t["attrs"][k]:
                            cand = t["attrs"][k][0]
                            break
                    prod = cand or "hypothetical protein"
                cds_id = f"cds-{tid}"
                cattrs_base = [
                    ("ID", cds_id),
                    ("Parent", rna_id),
                    ("gbkey", "CDS"),
                    ("product", prod),
                ]
                if dbx_t: cattrs_base.append(("Dbxref", dbx_t))
                for c in cds_sorted:
                    cattrs = list(cattrs_base)  # copy
                    emit(c["seqid"], args.source, "CDS", c["start"], c["end"], c["score"], c["strand"], c["frame"], cattrs)

    if out is not sys.stdout:
        out.close()

if __name__ == "__main__":
    main()
