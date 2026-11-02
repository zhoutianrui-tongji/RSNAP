#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

# 用户输入参数
user_args = {
    # 样本绝对路径（双末端）
    "samples": {  # 左                    #右
        "S1": ["./test/Con-1_1.fq.gz", "./test/Con-1_2.fq.gz"],
        "S2": ["./test/Con-2_1.fq.gz", "./test/Con-2_2.fq.gz"],
        "S3": ["./test/Sam-1_1.fq.gz", "./test/Sam-1_2.fq.gz"],
        "S4": ["./test/Sam-2_1.fq.gz", "./test/Sam-2_2.fq.gz"]
    },

    # 单末端文件示例：
    # "samples2": {
    #     "S1": ["./test/Con-1_1.fq.gz"],
    #     "S2": ["./test/Con-2_1.fq.gz"],
    #     "S3": ["./test/Sam-1_1.fq.gz"],
    #     "S4": ["./test/Sam-2_1.fq.gz"]
    # },

    # 样本分组信息
    "groups": {
        "group1": ["S1", "S2"],
        "group2": ["S3", "S4"]
    },

    # 测序类型（单双末端） paired或者single
    "layout": "paired",

    # 测序文件类型 fastq或者fasta
    "format": "fastq",

    # 结果目录（默认为./RSNAP_out/）
    "output_dir": "./RSNAP_out_new/",

    # 检索关键词，例如"(taxonomy_id:4751) AND (reviewed:true)"\n详细请查看手册：https://www.uniprot.org/help/query-fields
    "query": "(taxonomy_id:4751) AND (reviewed:true)"
}

# Quality_Control默认参数
qc_args = {
    "not_cut": False,  # 是否执行裁剪测序数据，默认为是(False表示不不切割，如过为True则为切割)
    "threads": "150",  # 线程数，默认为10
    "threshold": 0.01,  # 进行裁剪长度判断时的阈值，默认为0.01
    "software": "trimmomatic",  # 进行裁剪的程序，默认为trimmomatic
}

# Assemble_trinity默认参数
at_args = {
    "result_file": "assemblies.fa",  # 组装结果文件名
    "summary_file": "summary.txt",  # 组装摘要结果文件名
    "max_memory": "50",  # 最大内存占用数（G）（默认为6G）
    "cpu": "40",  # cup核数，默认为4
}

# Expression_evaluation_rsem默认参数
eer_args = {
    "threads": "100",  # 线程数，默认为8
}

# Differential_expression_gene_identification默认参数
degi_args = {
    "software": "edgeR",  # 差异分析R包,默认为'edgeR'。支持("edgeR", "DESeq2")
    "foldchange": "1",  # 差异分析log2foldchange筛选阈值，默认为1
    "padj": "0.05",  # 异分析矫正p值padj筛选阈值，默认为0.05
}

# Transcript_identification默认参数
ti_args = {
    "cpu": "40",  # CPU核数(默认为1)
    "evalue": 10 ** -5,  # evalue值(默认为10e-5)
}
# Feature_annotation_analysis默认参数
faa_args = {
    "all_out": "all_anno.tsv",  # 所有可鉴别转录本的各类注释信息结果文件名，默认为all_anno.tsv
    "deg_out": "deg_anno.tsv",  # 差异表达基因的各类注释信息结果文件名，默认为deg_anno.tsv
    "enrichment_out": "enrichment_out",  # 差异表达基因的各类注释信息的富集分析结果文件目录名,默认为enrichment_out（即在output目录下的enrichment_out子目录）
}
