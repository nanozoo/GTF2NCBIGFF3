# GTF → NCBI GFF3 Converter

A Python script to convert a **StringTie-like GTF** (or similar fungal genome annotation file)  
into an **NCBI GenBank-compatible GFF3** format.

This tool reconstructs the `gene → mRNA → exon` hierarchy and preserves/restructures relevant
functional annotations (Pfam, TIGR, InterPro, CDD, homology descriptions, GO terms, etc.)
into NCBI-friendly `Dbxref`, `Note`, `product`, `EC_number`, and related qualifiers — producing
output ready for GenBank/RefSeq submission.

---

## ✨ Features

- Converts `gene`, `transcript`, and `exon` entries to **NCBI-style** `gene`, `mRNA`, and `exon` features.
- Includes `CDS` features if present in your GTF (with correct phase).
- **Optionally infers CDS from exons** (`--infer-cds-from-exons`) for GTFs that lack CDS  
  *(enable only if your exons represent coding sequence only, i.e., no UTR exons).*
- Automatically constructs proper **GFF3 hierarchy** using `Parent=` and `ID=` attributes.
- Extracts and preserves useful metadata:
  - **Pfam / TIGRFAM / InterPro / CDD** → `Dbxref=...`
  - **UniProt accessions** (when detectable) → `Dbxref=UniProtKB:...`
  - **GO IDs** → `Dbxref=GO:...`
  - Pipeline annotations / HMMER matches / taxonomy / homologies, etc. → `Note=...`
- **Rescues functional CDS annotation wherever evidence exists**:
  - Uses `rec_name`, then `ref_tag`, then strong UniProt homologies (`sp|ACC|NAME`),  
    then HMM domain hints to populate meaningful `product=...`.
  - Adds `EC_number=...` when EC evidence is found.
  - Falls back to `hypothetical protein` only when no evidence is available.
- Properly escapes special characters (`;`, `,`, `=`) per GFF3 specification.
- Generates **systematic locus tags** like `FUNGI_00001`, `FUNGI_00002`, …
  - In **NCBI-strict mode**, `--locus-tag-prefix` is required and recommended for submission.
  - Use `--no-ncbi-strict` only for internal/testing output.
- Adds `protein_id` to CDS derived from locus tags (e.g., `FUNGI_00001.p1`).
- **Optional FASTA exports** from a genome FASTA:
  - CDS sequences (`--cds-fasta`)
  - Protein sequences (`--protein-fasta`)
  - FASTA headers include transcript + gene + locus_tag.
- If a genome FASTA is provided, emits `##sequence-region` pragmas automatically.
- Compatible with downstream NCBI pipelines such as **tbl2asn**, **tbl2asn_r10**, and GenBank submission tools.

---

## 🧩 Example

Input: `annotation.gtf` (from a genome annotation pipeline)  
Output: `annotation.ncbi.gff3` (NCBI-style)

```bash
# Basic usage (GFF3 to stdout)
python gtf_to_ncbi_gff3.py annotation.gtf > annotation.gff3
```

```bash
# Basic usage (write to file)
python gtf_to_ncbi_gff3.py annotation.gtf -o annotation.gff3
```

```bash
# NCBI-strict locus tags (recommended for submission)
python gtf_to_ncbi_gff3.py genome.gtf \
  --locus-tag-prefix FUNGI \
  -o genome.ncbi.gff3
```

```bash
# Customize the "source" field in the GFF3 (optional)
python gtf_to_ncbi_gff3.py genome.gtf \
  --locus-tag-prefix FUNGI \
  --source MyPipeline \
  -o genome.ncbi.gff3
```

```bash
# Infer CDS from exons (only if exons are coding-only)
python gtf_to_ncbi_gff3.py genome.gtf \
  --locus-tag-prefix FUNGI \
  --infer-cds-from-exons \
  -o genome.ncbi.gff3
```

```bash
# Provide genome FASTA to add ##sequence-region pragmas
# and export CDS + protein FASTA files
python gtf_to_ncbi_gff3.py genome.gtf \
  --locus-tag-prefix FUNGI \
  --genome-fasta genome.fna \
  --infer-cds-from-exons \
  --cds-fasta genome.cds.fna \
  --protein-fasta genome.prot.faa \
  -o genome.ncbi.gff3
```

```bash
# Use a different codon table (default 1 = Standard)
python gtf_to_ncbi_gff3.py genome.gtf \
  --locus-tag-prefix FUNGI \
  --genome-fasta genome.fna \
  --codon-table 1 \
  --protein-fasta genome.prot.faa \
  -o genome.ncbi.gff3
```

```bash
# Disable NCBI-strict mode (NOT recommended for submissions)
python gtf_to_ncbi_gff3.py genome.gtf \
  --no-ncbi-strict \
  -o genome.gff3
```

---

## 📦 Output Structure

Example snippet (protein-coding case):

```
##gff-version 3
##sequence-region contig_1 1 230218
contig_1  GTF2NCBI  gene  52524 56074 . - .  ID=gene-NANOZOOG1;Name=NANOZOOG1;locus_tag=FUNGI_00001;gbkey=Gene;Dbxref=...
contig_1  GTF2NCBI  mRNA  52524 56074 . - .  ID=rna-NANOZOOT1.1;Parent=gene-NANOZOOG1;gbkey=mRNA;transcript_id=NANOZOOT1.1;Dbxref=...;Note=...
contig_1  GTF2NCBI  exon  52524 53602 . - .  ID=exon-NANOZOOT1.1.1;Parent=rna-NANOZOOT1.1;gbkey=exon;exon_number=1
contig_1  GTF2NCBI  CDS   52580 56050 . - 0  ID=cds-NANOZOOT1.1;Parent=rna-NANOZOOT1.1;gbkey=CDS;product=Purine nucleoside permease;protein_id=FUNGI_00001.p1;EC_number=1.1.1.1;Dbxref=...
```

Each gene includes:

- `gene` → `mRNA` → `exon` (and optionally `CDS`)
- NCBI-compatible attributes:  
  `ID`, `Parent`, `locus_tag`, `Name`, `gbkey`, `Dbxref`, `Note`, `product`, `protein_id`, etc.

---

## ⚙️ Installation

Clone or download the script:

```bash
git clone https://github.com/nanozoo/GTF2NCBIGFF3.git
cd GTF2NCBIGFF3
```

Make executable (optional):

```bash
chmod +x gtf_to_ncbi_gff3.py
```

Run with Python 3.7+ (no external dependencies required):

```bash
python gtf_to_ncbi_gff3.py input.gtf -o output.gff3
```

---

## 🧠 How It Works

1. Parses each GTF line into structured data (regex captures `key "value"` attributes).
2. Groups features by `gene_id` and `transcript_id`.
3. Constructs hierarchical `gene → mRNA → exon/CDS` structure.
4. Extracts functional evidence:
   - Pfam / TIGRFAM / InterPro / CDD → `Dbxref`
   - UniProt accessions and GO IDs (when present) → `Dbxref`
   - Pipeline metadata / HMM matches / taxonomy / homology strings → `Note`
5. Infers meaningful CDS `product` where evidence exists (`rec_name`, `ref_tag`, strong UniProt hits, HMM domains).
6. Optionally infers CDS from exons and computes phase (`--infer-cds-from-exons`).
7. Encodes attributes safely for GFF3 output.
8. If genome FASTA is provided:
   - emits `##sequence-region` pragmas,
   - optionally exports CDS and protein FASTA files.

---

## 🧪 Requirements

- Python ≥ 3.7
- Works on Linux, macOS, and Windows (via WSL or Python terminal)
- Genome FASTA is **only required** if you want:
  - `##sequence-region` pragmas,
  - CDS/protein FASTA export,
  - or CDS inference sanity checks.

---

## 🧾 Notes for NCBI Submission

- For GenBank/RefSeq submission, always use a consistent `--locus-tag-prefix`.
- Validate your final GFF3 + genome FASTA with `tbl2asn` before submission.
- Only enable `--infer-cds-from-exons` if exons represent **coding sequence only**.
- NCBI may normalize overly specific product names; providing evidence-based names is still preferred.

---

## 📫 Feedback / Contributions

Pull requests and issues are welcome!  
If you encounter a special GTF dialect or NCBI validation edge case, please open an issue with a minimal example.
