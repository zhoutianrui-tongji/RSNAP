#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@File    :   Feature_annotation_analysis.py    
@Contact :   2454888366@qq.com

@Modify Time      @Author    @Version    @Desciption
------------      -------    --------    -----------
2023/3/31 上午10:00   skychou      1.0         None
'''
import os
import re
import logging as lg
import click
import gseapy as gp
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
matplotlib.use('Agg')
from gseapy import barplot, dotplot
import pickle

def get_identified_table(blast_res_file):
    with open(blast_res_file, "r") as f:
        blast_res_list = f.readlines()
    # 表头模式
    header_pat = re.compile("#(.*)BLAST(.*)")

    # 判断某行是否是表头
    def is_header(line):
        return re.findall(header_pat, line) != []

    # 获得表头在列表中索引
    header_index = []
    for n in range(len(blast_res_list)):
        line = blast_res_list[n]
        if is_header(line):
            header_index.append(n)
    # 获得每组blast比对结果
    every_blast = []
    for n in range(len(header_index) - 1):
        start = header_index[n] + 1
        end = header_index[n + 1]
        every_blast.append(blast_res_list[start:end])

    # 所有鉴别到的转录本
    resdict = {"query_id": [], "entry": []}
    for i in every_blast:
        if len(i) > 3:
            tempquery_id = i[4].split()[0]
            tempentry = i[4].split()[1]
            if tempquery_id not in resdict["query_id"]:
                resdict["query_id"].append(tempquery_id)
                resdict["entry"].append(tempentry)
    return pd.DataFrame(resdict)


def parse_blastlist(blast_list):
    return pd.read_csv(blast_list, sep='\t')


def get_tax_ids(tax_list):
    reslist = []
    for taxid in tax_list:
        match = re.search(r"\b(\d+)\s+\(kingdom\)", taxid)
        if match:
            reslist.append(match.group(1))
        else:
            reslist.append("")
    return reslist


def get_gene2anno(annotation_df, anno_type="Gene Ontology (GO)", outdir="./"):
    res_dict = {}
    tempdf = annotation_df.copy()
    df = tempdf[['Entry', anno_type]]
    df.to_csv(os.path.join(outdir, f"gene2anno_{anno_type}.tsv"), sep="\t", index=False)
    df = df.dropna()
    for i in range(len(df)):
        genes = str(df.iloc[i]["Entry"])
        genes = genes.split()
        annos = str(df.iloc[i][anno_type])
        annos = annos.split(";")
        for anno in annos:
            if anno not in res_dict:
                res_dict[anno] = []
            else:
                for gene in genes:
                    if gene not in res_dict[anno]:
                        res_dict[anno].append(gene)
    d = {k: v for k, v in res_dict.items() if v}
    return d


def get_out(annotation_file, df, output):
    annotation_df = pd.read_csv(annotation_file, sep='\t')
    annotation_df = annotation_df[annotation_df["Entry"].isin(df["entry"].tolist())]
    annotation_df.to_csv(output, sep="\t",index=False)
    return annotation_df


def format_filename(filename):
    # 将不规范的字符替换为下划线
    filename = re.sub(r'[^\w\s-]', '_', filename)
    # 将空格和连字符替换为下划线
    filename = filename.replace(' ', '_').replace('-', '_')
    # 去除开头和结尾的下划线
    filename = filename.strip('_')
    # 将所有字母转换为小写
    filename = filename.lower()
    # 去除多个下划线
    filename = re.sub(r'_+', '_', filename)
    return filename

#定义绘图和保存文件函数
def draw_write(pre,gene_list, gene_sets, anno_type, color,enrichment_out):
    enr = gp.enrich(gene_list=gene_list,  # or gene_list=glist
                    gene_sets=gene_sets,
                    background=background,  # "hsapiens_gene_ensembl",
                    verbose=True)
    if (enr.results['Adjusted P-value'] < 0.05).any():
        column = "Adjusted P-value"
    else:
        column = "P-value"
    ax = barplot(enr.results,
                 column=column,
                 group='Gene_set',  # set group, so you could do a multi-sample/library comparsion
                 title=f"{pre}_{anno_type}_barplot",
                 color=color
                 )
    ax.legend_.remove()
    ofname = os.path.join(enrichment_out, format_filename(f"{pre}_{anno_type}_barplot")+ ".pdf")
    fig = ax.figure  # 获取Figure对象
    fig.savefig(ofname, bbox_inches='tight', pad_inches=0)

    ax = dotplot(enr.results,
                 column=column,
                 x='Gene_set',  # set x axis, so you could do a multi-sample/library comparsion
                 title=f"{pre}_{anno_type}_dotplot",
                 show_ring=True,  # set to False to revmove outer ring
                 marker='o',
                 )
    ofname=os.path.join(enrichment_out, format_filename(f"{pre}_{anno_type}_dotplot") + ".pdf")
    ax.xaxis.set_visible(False)
    fig = ax.figure  # 获取Figure对象
    fig.savefig(ofname, bbox_inches='tight', pad_inches=0)

@click.command()
@click.option('--annotation_file', '-a', type=click.Path(), help='来自UniprotKB数据库的已知蛋白注释信息', prompt=True)
@click.option('--blast_outfile', '-b', type=click.Path(), help='BLAST比对鉴别结果文件', default="")
@click.option('--blast_list', '-l', type=click.Path(),
              help='BLAST比对鉴别结果文件解析后的基因列表文件（两列，列名为query_id和entry）', default="")
@click.option("--deg_list", "-d", help='差异表达基因列表文件（两列，列名为query_id和logFC）', default="")
@click.option('--output', '-o', type=click.Path(), help='结果输出文件(绝对路径)',
              default="./feature_annotation_analysis_output/")
@click.option('--all_out', '-ao', default="all_anno.tsv",
              help='所有可鉴别转录本的各类注释信息结果文件名，默认为all_anno.tsv')
@click.option('--deg_out', '-do', default="deg_anno.tsv",
              help='差异表达基因的各类注释信息结果文件名，默认为deg_anno.tsv')
@click.option('--enrichment_out', '-eo', default="enrichment_out",
              help='差异表达基因的各类注释信息的富集分析结果文件目录名,默认为enrichment_out（即在output目录下的enrichment_out子目录）')
@click.option("--config", "-c", type=click.Path(), default="",
              help="配置文件目录")
def main(annotation_file, blast_outfile, blast_list, deg_list, output, all_out, deg_out, enrichment_out, config):
    # blast_list = "./transcript_indentification_output/all_identified_transcripts.tsv"
    # deg_list = "./transcript_indentification_output/sigDEG.tsv"
    # annotation_file = "./taxonomy_id_4751_and_reviewed_true.summary"
    # output = './'
    # all_out = "all_anno.tsv"
    # deg_out = "deg_anno.tsv"

    # 获得config文件的内容
    global pipeline
    pipeline = False
    if config != "":
        # 获得config文件的内容
        global cfg
        if not os.path.exists(config):
            print(f"配置文件{config}不存在，请检查路径！")
            return
        try:
            import config as cfg
        except:
            sys.path.append(os.path.dirname(config))
            import config as cfg

        global config_prefix
        config_prefix = os.path.splitext(config)[0]
        pipeline = True

    # 设置日志文件
    # 记录器
    global logger
    logger = lg.getLogger("mylog")
    logger.setLevel(lg.DEBUG)

    # 处理器
    consleHandler = lg.StreamHandler()
    consleHandler.setLevel(lg.DEBUG)

    # out 目录检查并处理
    output_dir = output
    if output_dir == "./":
        output_dir = os.path.dirname(__file__)
    if output_dir != os.path.dirname(__file__):
        if output_dir.endswith("/"):
            output_dir = output_dir[:-1]
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
    if not os.path.exists(output_dir):
        raise Exception("out dir does not exist!输出路径不存在！")

    if not output_dir.endswith("/"):
        output_dir = output_dir + "/"
    output = output_dir

    # 没有给定handler日志级别，将用logger的级别
    fileHandler = lg.FileHandler(filename=f"{output_dir}feature_annotation_analysis.log", mode="w", encoding="UTF-8")
    fileHandler.setLevel(lg.INFO)

    # formatter格式
    # formatter = lg.Formatter("%(asctime)s|%(levelname)-8s|%(filename)10s:%(lineno)4s|%(message)s")
    formatter = lg.Formatter("%(asctime)s|%(filename)10s:%(lineno)4s|%(message)s")
    # 给处理器设置格式
    consleHandler.setFormatter(formatter)
    fileHandler.setFormatter(formatter)

    # 给记录器设置处理器
    logger.addHandler(consleHandler)
    logger.addHandler(fileHandler)


    # enrichment_out = "enrichment_out"
    if blast_outfile != "":
        blast_df = get_identified_table(blast_outfile)
    elif blast_list != "":
        blast_df = parse_blastlist(blast_list)
    else:
        logger.error(f"请给出必须参数，blast_outfile或者blast_list\n")
        raise Exception(f"请给出必须参数，blast_outfile或者blast_list")

    deg_df = pd.read_csv(deg_list, sep='\t')
    merged_df = pd.merge(deg_df, blast_df, on="query_id")

    annotation_df = pd.read_csv(annotation_file, sep='\t')

    # 获得差异表达基因的各类注释信息
    deg_out = os.path.join(output, deg_out)
    deg_out_df = get_out(annotation_file, merged_df, deg_out)

    # 获得所有可鉴别转录本的各类注释信息
    all_out = os.path.join(output, all_out)
    all_out_df = get_out(annotation_file, blast_df, all_out)

    # 获得背景信息：blast后的到的
    global background
    background = blast_df["entry"].tolist()
    background = list(set(background))

    # 富集分析与绘图
    enrichment_out = os.path.join(output, enrichment_out)
    if not os.path.exists(enrichment_out):
        os.makedirs(enrichment_out)

    # 上调基因：
    up_gene_list = merged_df[merged_df["logFC"] > 0]["entry"].tolist()
    up_gene_list = list(set(up_gene_list))
    # 下调基因：
    down_gene_list = merged_df[merged_df["logFC"] < 0]["entry"].tolist()
    down_gene_list = list(set(down_gene_list))

    #先进行GO富集分析
    anno_type = "Gene Ontology (GO)"
    logger.info(f"开始进行{anno_type}富集分析\n")
    gene_sets = get_gene2anno(annotation_df, anno_type, outdir=output)
    #上调：
    draw_write("up",up_gene_list,gene_sets,anno_type,"red",enrichment_out)
    #下调
    draw_write("down", down_gene_list, gene_sets, anno_type, "green", enrichment_out)
    logger.info(f"{anno_type}富集分析完成！\n")
    #进行kegg通路富集分析
    anno_type = "Pathway"
    logger.info(f"开始进行{anno_type}富集分析\n")
    annotation_df[anno_type] = annotation_df[anno_type].str.replace('PATHWAY: ', '')
    gene_sets = get_gene2anno(annotation_df, anno_type, outdir=output)
    # 上调：
    draw_write("up", up_gene_list, gene_sets, anno_type, "red", enrichment_out)
    # 下调
    draw_write("down", down_gene_list, gene_sets, anno_type, "green", enrichment_out)
    logger.info(f"{anno_type}富集分析完成！\n")
    # 读取原始 pkl 文件中的列表数据
    if pipeline:
        config_pkl = f'{config_prefix}.pkl'
        with open(config_pkl, "rb") as f:
            data = pickle.load(f)

        # 在原始列表中添加新的字符串元素
        enrichment_out_file = os.listdir(enrichment_out)    # 获取out下的所有文件和文件夹
        res_dict = {"all_out":all_out,"deg_out":deg_out,"enrichment_out":enrichment_out_file}
        data["annotation_analysis_files"] = res_dict

        # 将原始列表和新的字符串一起写入新的 pkl 文件中
        with open(config_pkl, "wb") as f:
            pickle.dump(data, f)
    logger.info(f"特征注释分析完成！\n")


if __name__ == '__main__':
    main()
