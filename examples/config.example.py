"""Example RSNAP configuration.

Copy this file, replace the sample paths and group labels, then pass it to
script/RSNAP.py with --config.
"""

user_args = {
    "samples": {
        "control_1": ["/path/to/control_1_R1.fastq.gz", "/path/to/control_1_R2.fastq.gz"],
        "control_2": ["/path/to/control_2_R1.fastq.gz", "/path/to/control_2_R2.fastq.gz"],
        "treated_1": ["/path/to/treated_1_R1.fastq.gz", "/path/to/treated_1_R2.fastq.gz"],
        "treated_2": ["/path/to/treated_2_R1.fastq.gz", "/path/to/treated_2_R2.fastq.gz"],
    },
    "groups": {
        "control": ["control_1", "control_2"],
        "treated": ["treated_1", "treated_2"],
    },
    "layout": "paired",
    "format": "fastq",
    "output_dir": "./RSNAP_out",
    "query": "(taxonomy_id:4751) AND (reviewed:true)",
    "addition_fields": "",
    "compare": [["control", "treated"]],
}

qc_args = {
    "not_cut": False,
    "threads": "8",
    "threshold": 0.01,
    "software": "trimmomatic",
    "script": "./Quality_Control.py",
}

at_args = {
    "result_file": "assemblies.fa",
    "summary_file": "summary.txt",
    "max_memory": "16",
    "cpu": "8",
    "script": "./Assemble_trinity.py",
}

eer_args = {"threads": "8", "script": "./Expression_evaluation_rsem.py"}
degi_args = {
    "software": "edgeR",
    "foldchange": "1",
    "padj": "0.05",
    "script": "./Differential_expression_gene_identification.py",
}
ti_args = {
    "cpu": "8",
    "evalue": 1e-5,
    "script": "./Transcript_identification.py",
}
faa_args = {
    "all_out": "all_anno.tsv",
    "deg_out": "deg_anno.tsv",
    "enrichment_out": "enrichment_out",
    "script": "./Feature_annotation_analysis.py",
}
