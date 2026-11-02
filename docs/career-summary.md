# RSNAP career summary

## Chinese resume version

**RSNAP：RNA-seq 无参考分析自动化管道｜Python / R / Bioinformatics**

- 独立开发 RNA-seq 无参考分析流程，整合 FastQC、Trimmomatic、Trinity、RSEM、edgeR/DESeq2、BLASTX、UniProt 和 GSEApy。
- 实现从原始 FASTQ 质控、读段裁剪、转录组组装、表达定量、差异表达分析到功能注释和 GO/通路富集分析的自动化执行。
- 设计配置驱动的样本分组和参数管理机制，记录阶段性输出路径，并生成统一 HTML 分析报告。
- 在真实 RNA-seq 数据上完成流程验证，示例运行获得 24,174 条组装转录本和 1,143 个显著差异表达基因。

## English resume version

**RSNAP: Automated de novo RNA-seq analysis pipeline | Python / R / Bioinformatics**

- Developed a configurable de novo RNA-seq workflow integrating FastQC, Trimmomatic, Trinity, RSEM, edgeR/DESeq2, BLASTX, UniProt, and GSEApy.
- Automated read QC, trimming, transcriptome assembly, expression quantification, differential expression analysis, transcript identification, and enrichment analysis.
- Implemented configuration-driven sample grouping, intermediate-result tracking, and consolidated HTML report generation.
- Validated the workflow on real RNA-seq data, producing 24,174 assembled transcripts and 1,143 significant differentially expressed genes in an archived run.

## Interview explanation

The project demonstrates end-to-end workflow engineering in a scientific computing setting. The key engineering work was translating a multi-tool analysis protocol into a reproducible pipeline with explicit inputs, outputs, configuration, logging, and report generation. The individual tools provide the underlying statistical and sequence-analysis methods; RSNAP connects them into one usable workflow.
