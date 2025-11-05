# GTF → NCBI GFF3 Converter

A Python script to convert a **StringTie-like GTF** (or similar (fungal) genome annotation file)  
into an **NCBI GenBank-compatible GFF3** format.

This tool reconstructs the `gene → mRNA → exon` hierarchy and preserves relevant
functional annotations (Pfam, TIGR, InterPro, HMMER descriptions, etc.)
in `Dbxref` and `Note` fields—producing output ready for GenBank/RefSeq submission.

## ✨ Features

- Converts `gene`, `transcript`, and `exon` entries to **NCBI-style** `gene`, `mRNA`, and `exon` features.
- Optionally includes `CDS` features if present in your GTF (with correct phase).
- Automatically constructs proper **GFF3 hierarchy** using `Parent=` and `ID=` attributes.
- Extracts and preserves useful metadata:
  - **Pfam / TIGRFAM / InterPro / CDD** → `Dbxref=...`
  - **HMMER & desc fields** → `Note=...`
- Properly escapes special characters (`;`, `,`, `=`) per GFF3 specification.
- Optionally generates **standardized locus tags** like `FUNGI_00001`, `FUNGI_00002`, …
- Compatible with downstream NCBI pipelines such as **tbl2asn**, **tbl2asn_r10**, and **GenBank submission tools**.

## 🧩 Example

Input: `annotation.gtf` (from a genome annotation pipeline)  
Output: `annotation.gff3` (NCBI-style)

```bash
# Basic usage
python gtf_to_ncbi_gff3.py annotation.gtf -o annotation.gff3
````

```bash
# Add a locus tag prefix (auto-numbered genes)
python gtf_to_ncbi_gff3.py genome.gtf -o genome.gff3 --locus-tag-prefix FUNGI
```

```bash
# Customize the "source" field in the GFF3 (optional)
python gtf_to_ncbi_gff3.py genome.gtf -o genome.gff3 --source MyPipeline
```

---

## 📦 Output Structure

Example snippet:

```
##gff-version 3
contig_1_segment0_pilon_pilon	GTF2NCBI	gene	52524	56074	1000	-	.	ID=gene-NANOZOOG1;Name=NANOZOOG1;locus_tag=NANOZOOG1;gbkey=Gene
contig_1_segment0_pilon_pilon	GTF2NCBI	mRNA	52524	56074	1000	-	.	ID=rna-NANOZOOT1.1;Parent=gene-NANOZOOG1;gbkey=mRNA;Dbxref=Pfam:PF00533.22,Pfam:PF03031.14;Note=HMMER_1_desc=BRCT...
contig_1_segment0_pilon_pilon	GTF2NCBI	exon	52524	53602	1000	-	.	ID=exon-NANOZOOT1.1.1;Parent=rna-NANOZOOT1.1;gbkey=mRNA;exon_number=1
```

Each gene includes:

* `gene` → `mRNA` → `exon` (and optionally `CDS`)
* NCBI-compatible attributes:
  `ID`, `Parent`, `locus_tag`, `Name`, `gbkey`, `Dbxref`, `Note`, `product`, etc.

---

## ⚙️ Installation

Clone or download the script:

```bash
git clone https://github.com/<yourusername>/gtf-to-ncbi-gff3.git
cd gtf-to-ncbi-gff3
```

Make executable:

```bash
chmod +x gtf_to_ncbi_gff3.py
```

Run with Python 3.7+ (no external dependencies required):

```bash
python gtf_to_ncbi_gff3.py input.gtf -o output.gff3
```

---

## 🧠 How It Works

1. Parses each GTF line into structured data (using regex to capture attributes).
2. Groups features by `gene_id` and `transcript_id`.
3. Constructs hierarchical `gene → mRNA → exon/CDS` structure.
4. Extracts relevant description fields:
   * Pfam / TIGRFAM / InterPro → `Dbxref`
   * HMMER / description / match fields → `Note`
5. Encodes attributes safely for GFF3 output.

---

## 🧪 Requirements

* Python ≥ 3.7
* Works on Linux, macOS, and Windows (via WSL or Python terminal)

