#!/usr/bin/env python3
"""
Convert fungal GTF annotations to NCBI-compliant GFF3.

Features:
- Builds gene -> mRNA -> exon hierarchy.
- Preserves CDS if present in GTF.
- Optionally infers CDS from exons (--infer-cds-from-exons).
- Enforces systematic locus_tag with registered prefix (NCBI strict mode).
- Adds protein_id to CDS derived from locus_tag.
- Rescues functional annotations from rich GTF attributes:
    * product from rec_name/ref_tag/strong homology/hmm_matches
    * EC numbers into EC_number
    * GO IDs + UniProt/Pfam/InterPro/CDD into Dbxref
    * other metadata into Note
- Optionally exports CDS and protein FASTA from a genome FASTA.
- If genome FASTA is provided, emits ##sequence-region pragmas.

Usage:
  python gtf_to_ncbi_gff3.py input.gtf -o out.gff3 --locus-tag-prefix ABC
  python gtf_to_ncbi_gff3.py input.gtf --locus-tag-prefix ABC \
      --genome-fasta genome.fna --infer-cds-from-exons \
      --cds-fasta out.cds.fna --protein-fasta out.prot.faa
"""

import sys, re, argparse, collections
from urllib.parse import quote

# ------------------ helpers ------------------

def parse_gtf_attributes(attr_field: str):
    d = collections.defaultdict(list)
    for m in re.finditer(r'(\S+)\s+"([^"]*)"', attr_field):
        d[m.group(1)].append(m.group(2))
    d["_raw"] = [attr_field]
    return d

def first(d, *keys, default=None):
    for k in keys:
        if k in d and d[k]:
            return d[k][0]
    return default

def gff3_escape(v: str) -> str:
    # GFF3 percent-encoding. Keep common safe ASCII.
    return quote(v, safe=":+._-")

def join_attrs(pairs):
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

def collect_dbxrefs(raw_attrs: str):
    """Pull common domain/database IDs out of raw attribute text."""
    xrefs = set()
    for acc in re.findall(r"(Pfam:PF\d+(?:\.\d+)?)", raw_attrs):
        xrefs.add(acc)
    for acc in re.findall(r"(TIGRFAMs?:TIGR\d+)", raw_attrs):
        xrefs.add(acc.replace("TIGRFAMs:", "TIGRFAM:"))
    for acc in re.findall(r"(CDD:cd\d+)", raw_attrs, flags=re.IGNORECASE):
        xrefs.add(acc.upper().replace("CDD:CD", "CDD:cd"))
    for acc in re.findall(r"(InterPro:IPR\d+)", raw_attrs):
        xrefs.add(acc)
    for acc in re.findall(r"Dbxref=['\"]([^'\"]+)['\"]", raw_attrs):
        xrefs.add(acc)
    return sorted(xrefs)

def summarize_notes(attrs):
    """Put non-standard but useful fields into Note."""
    note_bits = []
    for key in sorted(attrs.keys()):
        if key.startswith("_"):
            continue
        low = key.lower()
        if any(tok in low for tok in ["desc", "hmm", "match", "note",
                                      "homolog", "annotation_tag", "taxonomy",
                                      "source"]):
            vals = list(collections.OrderedDict.fromkeys(attrs[key]))
            s = f"{key}=" + ",".join(vals[:3])
            if len(vals) > 3:
                s += f",(+{len(vals)-3} more)"
            note_bits.append(s)
    n = "; ".join(note_bits)
    return n[:900]

# ---------- functional rescue ----------

def clean_uniprot_name(s):
    return s.replace("_", " ").strip()

def extract_product_from_attrs(attrs):
    """
    Conservative product extraction:
    1) rec_name Full=
    2) ref_tag
    3) strong UniProt homology sp|ACC|NAME
    4) hmm_matches domain-based (skip NA)
    else None
    """
    # 1) rec_name Full=...
    rec = first(attrs, "rec_name")
    if rec:
        m = re.search(r"Full=([^,{]+)", rec)
        if m:
            prod = re.sub(r"\{[^}]*\}", "", m.group(1)).strip()
            if prod:
                return clean_uniprot_name(prod)

    # 2) ref_tag
    ref_tag = first(attrs, "ref_tag")
    if ref_tag:
        return clean_uniprot_name(ref_tag)

    # 3) strong UniProt style homology sp|ACC|NAME_SPECIES
    hom = first(attrs, "homologies")
    if hom:
        m = re.search(r"\b(?:sp|tr)\|([A-Z0-9]+)\|([^:,\s]+)", hom)
        if m:
            return clean_uniprot_name(m.group(2))

    # 4) hmm_matches -> domain-based name
    hmm = first(attrs, "hmm_matches")
    if hmm and hmm.upper() not in {"NA", "N/A", ".", ""}:
        top_dom = hmm.split(":")[0].split(",")[0]
        return f"protein with {top_dom} domain"

    return None

def extract_ec_numbers(attrs):
    rec = first(attrs, "rec_name")
    if not rec:
        return []
    return re.findall(r"EC=([0-9]+\.[0-9]+\.[0-9]+\.[0-9]+)", rec)

def extract_go_ids(attrs):
    gos = []
    go_raw = first(attrs, "GO")
    if go_raw:
        gos.extend(re.findall(r"(GO:\d{7})", go_raw))
    return sorted(set(gos))

def extract_uniprot_accessions(attrs):
    """
    Extract UniProt accessions from homologies.
    - strong sp|ACC| / tr|ACC|
    - weak ACC_SPECIES (Dbxref only, not for product)
    """
    accs = []
    hom = first(attrs, "homologies")
    if not hom:
        return []
    accs += re.findall(r"\b(?:sp|tr)\|([A-Z0-9]+)\|", hom)
    accs += re.findall(r"\b([A-NR-Z0-9]{6,10})_[A-Z0-9]{3,5}\b", hom)
    return sorted(set(accs))

# ---------- FASTA / translation helpers ----------

def load_fasta(path):
    seqs = {}
    name = None
    buf = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf).upper()
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if name:
            seqs[name] = "".join(buf).upper()
    return seqs

_comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def revcomp(s): return s.translate(_comp)[::-1]

CODON_TABLES = {
    1: {  # Standard
        **{k:v for k,v in {
            "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
            "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
            "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
            "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
            "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
            "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
            "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
            "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
            "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
            "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
            "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
            "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
            "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
            "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
            "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
            "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
        }.items()}
    }
}

def translate(seq, table=1, to_stop=True):
    codons = CODON_TABLES.get(table, CODON_TABLES[1])
    prot = []
    for i in range(0, len(seq) - 2, 3):
        aa = codons.get(seq[i:i+3], "X")
        if aa == "*" and to_stop:
            break
        prot.append(aa)
    return "".join(prot)

def compute_cds_phases(segments):
    """Compute phase for inferred CDS segments (transcription order)."""
    phases = []
    offset = 0
    for (s,e) in segments:
        length = e - s + 1
        phase = offset % 3
        phases.append(str(phase))
        offset += length
    return phases

# ------------------ main ------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gtf", help="Input GTF")
    ap.add_argument("-o", "--out", help="Output GFF3 file (default: stdout)")
    ap.add_argument("--source", default="GTF2NCBI", help="Source field to use in output")
    ap.add_argument("--locus-tag-prefix", dest="lt_prefix", default=None,
                    help="NCBI locus_tag prefix (e.g. ABC). Required in strict mode.")
    ap.add_argument("--no-ncbi-strict", action="store_true",
                    help="Allow missing locus-tag prefix (NOT recommended for NCBI).")
    ap.add_argument("--default-product", default="hypothetical protein",
                    help="Fallback product for CDS if no evidence is found")
    ap.add_argument("--infer-cds-from-exons", action="store_true",
                    help="If no CDS in GTF, treat exons as CDS and compute phase.")
    ap.add_argument("--genome-fasta", help="Genome FASTA for CDS/protein export and sequence-region pragmas")
    ap.add_argument("--cds-fasta", help="Output CDS FASTA")
    ap.add_argument("--protein-fasta", help="Output protein FASTA")
    ap.add_argument("--codon-table", type=int, default=1,
                    help="NCBI codon table id (default 1)")
    args = ap.parse_args()

    if not args.no_ncbi_strict and not args.lt_prefix:
        ap.error("--locus-tag-prefix is required for NCBI-compliant output.")

    genome = None
    if args.genome_fasta:
        genome = load_fasta(args.genome_fasta)

    out = sys.stdout if not args.out else open(args.out, "w", encoding="utf-8")
    cds_out = open(args.cds_fasta, "w", encoding="utf-8") if args.cds_fasta else None
    prot_out = open(args.protein_fasta, "w", encoding="utf-8") if args.protein_fasta else None

    # Header
    print("##gff-version 3", file=out)

    # sequence-region pragmas if genome FASTA known
    if genome:
        for seqid in sorted(genome.keys()):
            print(f"##sequence-region {seqid} 1 {len(genome[seqid])}", file=out)

    genes = {}
    transcripts = collections.defaultdict(lambda: {"exons": [], "cds": []})

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

            if not gid:
                continue

            feature_low = feature.lower()
            if feature_low == "gene":
                genes[gid] = {
                    "seqid": seqid, "start": int(start), "end": int(end),
                    "score": score if score != "." else ".",
                    "strand": strand, "attrs": attrs,
                }
            elif feature_low in ("transcript", "mrna"):
                transcripts[(gid, tid)].update({
                    "seqid": seqid, "start": int(start), "end": int(end),
                    "score": score if score != "." else ".",
                    "strand": strand, "attrs": attrs,
                })
            elif feature_low == "exon":
                transcripts[(gid, tid)]["exons"].append({
                    "seqid": seqid, "start": int(start), "end": int(end),
                    "score": score if score != "." else ".",
                    "strand": strand, "attrs": attrs,
                    "exon_number": first(attrs, "exon_number", "exonNumber")
                })
            elif feature_low == "cds":
                phase = frame if frame in ("0","1","2") else "."
                transcripts[(gid, tid)]["cds"].append({
                    "seqid": seqid, "start": int(start), "end": int(end),
                    "score": score if score != "." else ".",
                    "strand": strand, "frame": phase, "attrs": attrs
                })

    # Group transcripts per gene
    transcripts_by_gene = collections.defaultdict(list)
    for (gid, tid), t in transcripts.items():
        if tid is None:
            continue
        transcripts_by_gene[gid].append((tid, t))

    # ID uniqueness check
    used_ids = set()
    def require_unique(id_):
        if id_ in used_ids:
            raise ValueError(f"Duplicate ID detected: {id_}")
        used_ids.add(id_)

    def emit(seqid, type_, start, end, score, strand, phase, attrs_pairs):
        print("\t".join([
            str(seqid), args.source, type_, str(start), str(end),
            score, strand, phase, join_attrs(attrs_pairs)
        ]), file=out)

    counter = 0

    for gid, g in genes.items():
        counter += 1
        seqid = g["seqid"]; strand = g["strand"]

        gene_id = f"gene-{gid}"
        require_unique(gene_id)

        locus_tag = f"{args.lt_prefix}_{counter:05d}" if args.lt_prefix else gid
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

        emit(seqid, "gene", g["start"], g["end"], g["score"], strand, ".", gene_attrs)

        # transcripts for this gene
        for iso_idx, (tid, t) in enumerate(sorted(transcripts_by_gene.get(gid, [])), start=1):

            rna_id = f"rna-{tid}"
            require_unique(rna_id)

            # If transcript bounds missing, infer from exons
            if not t.get("start") or not t.get("end"):
                if t["exons"]:
                    starts = [e["start"] for e in t["exons"]]
                    ends = [e["end"] for e in t["exons"]]
                    t["start"], t["end"] = min(starts), max(ends)

            traw = t.get("attrs", {}).get("_raw", [""])[0]
            dbx_t = sorted(set(dbx) | set(collect_dbxrefs(traw)))
            note_t = summarize_notes(t.get("attrs", {}))

            mrna_attrs = [
                ("ID", rna_id),
                ("Parent", gene_id),
                ("Name", tid),
                ("gbkey", "mRNA"),
                ("transcript_id", tid),
            ]
            if dbx_t: mrna_attrs.append(("Dbxref", dbx_t))
            if note_t: mrna_attrs.append(("Note", note_t))

            emit(t["seqid"], "mRNA", t["start"], t["end"], t["score"], t["strand"], ".", mrna_attrs)

            # Exons
            exons_sorted = sorted(
                t["exons"], key=lambda e: (e["start"], e["end"]),
                reverse=(t["strand"]=="-")
            )
            for i, e in enumerate(exons_sorted, start=1):
                ex_num = e.get("exon_number") or str(i)
                exon_id = f"exon-{tid}.{ex_num}"
                require_unique(exon_id)

                eattrs = [
                    ("ID", exon_id),
                    ("Parent", rna_id),
                    ("gbkey", "exon"),
                    ("exon_number", ex_num),
                ]
                emit(e["seqid"], "exon", e["start"], e["end"], e["score"], e["strand"], ".", eattrs)

            # CDS handling
            cds_segments = t["cds"]
            inferred = False

            if not cds_segments and args.infer_cds_from_exons and exons_sorted:
                cds_segments = [
                    {"seqid": e["seqid"], "start": e["start"], "end": e["end"],
                     "score": e["score"], "strand": e["strand"], "frame": "."}
                    for e in exons_sorted
                ]
                inferred = True

            if cds_segments:
                cds_sorted = sorted(
                    cds_segments, key=lambda c: (c["start"], c["end"]),
                    reverse=(t["strand"]=="-")
                )

                if inferred:
                    seg_coords = [(c["start"], c["end"]) for c in cds_sorted]
                    phases = compute_cds_phases(seg_coords)
                    for c, ph in zip(cds_sorted, phases):
                        c["frame"] = ph

                prod = extract_product_from_attrs(t.get("attrs", {}))
                if not prod:
                    prod = args.default_product

                ecs = extract_ec_numbers(t.get("attrs", {}))
                gos = extract_go_ids(t.get("attrs", {}))
                uniprot_accs = extract_uniprot_accessions(t.get("attrs", {}))

                dbx_extra = []
                for a in uniprot_accs:
                    dbx_extra.append(f"UniProtKB:{a}")
                for g_id in gos:
                    dbx_extra.append(g_id)

                dbx_all = sorted(set(dbx_t) | set(dbx_extra)) if (dbx_t or dbx_extra) else []

                cds_id = f"cds-{tid}"
                require_unique(cds_id)

                protein_id = f"{locus_tag}.p{iso_idx}"

                cattrs_base = [
                    ("ID", cds_id),
                    ("Parent", rna_id),
                    ("gbkey", "CDS"),
                    ("product", prod),
                    ("protein_id", protein_id),
                    ("transcript_id", tid),
                ]
                if ecs:
                    cattrs_base.append(("EC_number", ecs))
                if dbx_all:
                    cattrs_base.append(("Dbxref", dbx_all))

                for c in cds_sorted:
                    emit(c["seqid"], "CDS", c["start"], c["end"],
                         c.get("score", "."), c.get("strand", strand),
                         c.get("frame", "."), list(cattrs_base))

                # FASTA exports
                if genome and (cds_out or prot_out):
                    try:
                        seq = "".join(genome[c["seqid"]][c["start"]-1:c["end"]] for c in cds_sorted)
                    except KeyError as e:
                        raise KeyError(f"Contig {e} not found in genome FASTA.") from e

                    if t["strand"] == "-":
                        seq = revcomp(seq)

                    if cds_out:
                        cds_out.write(f">{tid}|{gid}|{locus_tag}\n")
                        for j in range(0, len(seq), 60):
                            cds_out.write(seq[j:j+60] + "\n")

                    if prot_out:
                        prot = translate(seq, table=args.codon_table, to_stop=True)
                        prot_out.write(f">{protein_id}|{tid}|{gid}|{locus_tag}\n")
                        for j in range(0, len(prot), 60):
                            prot_out.write(prot[j:j+60] + "\n")

    if out is not sys.stdout:
        out.close()
    if cds_out:
        cds_out.close()
    if prot_out:
        prot_out.close()

if __name__ == "__main__":
    main()
