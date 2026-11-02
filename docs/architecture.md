# RSNAP architecture

RSNAP is organized as a Python orchestrator around established command-line and R-based bioinformatics tools. Each stage writes files into a dedicated output directory and records paths for downstream stages.

```mermaid
flowchart LR
    A[FASTQ reads] --> B[Quality_Control.py]
    B --> B1[FastQC metrics]
    B --> B2[Trimmed reads]
    B2 --> C[Assemble_trinity.py]
    C --> C1[assemblies.fa]
    C --> C2[summary.txt]
    C1 --> D[Expression_evaluation_rsem.py]
    B2 --> D
    D --> D1[genes.expected_count.matrix]
    D1 --> E[Differential_expression_gene_identification.py]
    E --> E1[DEG tables and plots]
    C1 --> F[Transcript_identification.py]
    F --> F1[BLASTX and UniProt results]
    E1 --> G[Feature_annotation_analysis.py]
    F1 --> G
    G --> G1[GO/pathway enrichment]
    B1 --> H[HTML report]
    C2 --> H
    D1 --> H
    E1 --> H
    F1 --> H
    G1 --> H
```

## Control flow

`script/RSNAP.py` reads the configuration, creates output directories, invokes each stage, checks for existing outputs, and assembles the final Jinja2 report. The stage scripts can also be run independently through their Click command-line interfaces.

## Data flow

The workflow keeps the biological data units explicit: reads are used for quality control and assembly; assembled transcripts are used for quantification and BLASTX identification; count matrices are used for differential analysis; annotated transcript IDs are used for enrichment analysis.

## External tools

FastQC and Trimmomatic handle read quality, Trinity performs de novo assembly, RSEM estimates expression, edgeR or DESeq2 performs differential analysis, BLAST+ and UniProt provide protein-level identification, and GSEApy produces enrichment results and plots.
