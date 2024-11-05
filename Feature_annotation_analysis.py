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
import matplotlib.cm as cm
import matplotlib.colors as Colors
import numpy as np

matplotlib.use('Agg')  # 或者使用其他交互式后端
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
    df = tempdf[['query_id', anno_type]]
    df.to_csv(os.path.join(outdir, f"gene2anno_{anno_type}.tsv"), sep="\t", index=False)
    df = df.dropna()
    for i in range(len(df)):
        genes = str(df.iloc[i]["query_id"])
        genes = genes.split()
        annos = str(df.iloc[i][anno_type])
        annos = annos.split(";")
        for anno in annos:
            anno = re.sub(r'\s*\[GO:\d+\]', '', anno)
            anno = anno.strip()
            if anno not in res_dict:
                res_dict[anno] = []
            else:
                for gene in genes:
                    if gene not in res_dict[anno]:
                        res_dict[anno].append(gene)
    d = {k: v for k, v in res_dict.items() if v}
    return d


def get_go_type(annotation_df):
    tempdf = annotation_df.copy()
    df = tempdf[['Gene Ontology (biological process)', "Gene Ontology (cellular component)",
                 "Gene Ontology (molecular function)"]]
    bp = []
    cc = []
    mf = []
    for i in range(len(df)):
        annos = df.iloc[i]['Gene Ontology (biological process)']
        if not pd.isnull(annos):
            annos = annos.split(';')
            for anno in annos:
                anno = re.sub(r'\s*\[GO:\d+\]', '', anno)
                anno = anno.strip()
                if anno not in bp:
                    bp.append(anno)
        annos = df.iloc[i]["Gene Ontology (cellular component)"]
        if not pd.isnull(annos):
            annos = annos.split(';')
            for anno in annos:
                anno = re.sub(r'\s*\[GO:\d+\]', '', anno)
                anno = anno.strip()
                if anno not in cc:
                    cc.append(anno)
        annos = df.iloc[i]["Gene Ontology (molecular function)"]
        if not pd.isnull(annos):
            annos = annos.split(';')
            for anno in annos:
                anno = re.sub(r'\s*\[GO:\d+\]', '', anno)
                anno = anno.strip()
                if anno not in mf:
                    mf.append(anno)
    Term = bp + cc + mf
    go_type = ["BP"] * len(bp) + ["CC"] * len(cc) + ["MF"] * len(mf)
    data = {'Term': Term, 'go_type': go_type}
    df = pd.DataFrame(data)
    return df


def get_out(annotation_df, df, output):
    annotation_df = annotation_df[annotation_df["query_id"].isin(df["query_id"].tolist())]
    annotation_df.to_csv(output, sep="\t", index=False)
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


# 定义绘图和保存文件函数
def draw_write(pre, gene_list, gene_sets, anno_type, n_count, color_bar, enrichment_out):
    enr = gp.enrich(gene_list=gene_list,  # or gene_list=glist
                    gene_sets=gene_sets,
                    background=background,  # "hsapiens_gene_ensembl",
                    verbose=True)
    count = (enr.results['Adjusted P-value'] < 0.05).sum()
    if (enr.results['Adjusted P-value'] < 0.05).any() and count > 5:
        column = "Adjusted P-value"
    else:
        column = "P-value"

    df = enr.results
    df['count'] = df['Genes'].str.count(';') + 1
    df.sort_values(by='count', inplace=True, ascending=False)
    df['Term'] = df['Term'].apply(lambda x: x.strip())
    df = pd.merge(df, go_type_df, on='Term', how='inner')

    # df = enr.results
    # df['count'] = df['Genes'].str.count(';') + 1
    # df.sort_values(by='count', inplace=True, ascending=False)

    # 创建自定义颜色映射
    cmap = cm.get_cmap(color_bar)  # 使用红绿渐变色阶
    norm = Colors.Normalize(vmin=df[column].min(),
                            vmax=df[column].max())
    colors = [cmap(norm(value)) for value in df[column][0:n_count]]

    # 绘制柱状图
    fig, ax = plt.subplots()
    bars = ax.barh(df['Term'][0:n_count], df['count'][0:n_count], color=colors)

    # 创建颜色条
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # 必须设置一个空数组

    cbar = plt.colorbar(sm)
    cbar.set_label(column)

    # 设置图形属性
    ax.set_xlabel('Count')
    ax.set_ylabel('Term')
    ax.set_title(f"{pre}_{anno_type}_barplot")

    ofname = os.path.join(enrichment_out, format_filename(f"{pre}_{anno_type}_barplot") + ".pdf")
    plt.savefig(ofname, bbox_inches='tight', pad_inches=0.2)
    plt.clf()

    # 计算气泡图的大小和位置
    x = (df["Odds Ratio"][0:n_count])
    y = df['Term'][0:n_count]
    sizes = (df["count"][0:n_count] * 20)

    # 绘制气泡图
    fig, ax = plt.subplots()
    scatter = ax.scatter(x, y, s=sizes, c=colors)

    # 创建颜色条
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # 必须设置一个空数组

    cbar = plt.colorbar(sm)
    cbar.set_label(column)

    # 设置图形属性
    ax.set_xlabel('Odds Ratio')
    ax.set_title(f"{pre}_{anno_type}_dotplot")

    # 增加表示气泡大小的图例
    handles, labels = scatter.legend_elements(prop="sizes", alpha=0.1, func=lambda x: x / 20)
    legend2 = ax.legend(handles, labels, loc="upper right", title="Count", bbox_to_anchor=(1.4, 1), borderaxespad=0)

    ofname = os.path.join(enrichment_out, format_filename(f"{pre}_{anno_type}_dotplot") + ".pdf")
    plt.savefig(ofname, bbox_inches='tight', pad_inches=0.2)
    plt.clf()
    # ofname=os.path.join(enrichment_out, format_filename(f"{pre}_{anno_type}_dotplot") + ".pdf")
    # fig = ax.figure  # 获取Figure对象
    # fig.savefig(ofname, bbox_inches='tight', pad_inches=0)


def draw_write_go(pre, gene_list, gene_sets, anno_type, n_count, color_bar, enrichment_out):
    enr = gp.enrich(gene_list=gene_list,  # or gene_list=glist
                    gene_sets=gene_sets,
                    background=background,  # "hsapiens_gene_ensembl",
                    verbose=True)
    count = (enr.results['Adjusted P-value'] < 0.05).sum()
    if (enr.results['Adjusted P-value'] < 0.05).any() and count > 5:
        column = "Adjusted P-value"
    else:
        column = "P-value"

    df = enr.results
    df['count'] = df['Genes'].str.count(';') + 1

    # 自定义聚合函数，用于将字符串相加
    def concat_strings(x):
        return ';'.join(x)

    grouped_df = df.groupby('Term').agg({
        'Gene_set': 'first',
        'Overlap': 'sum',
        'P-value': 'mean',
        'Adjusted P-value': 'mean',
        'Odds Ratio': 'mean',
        'Genes': concat_strings,
        'count': 'sum',
    }).reset_index()
    df = grouped_df
    df['Term'] = df['Term'].apply(lambda x: x.strip())
    df = pd.merge(df, go_type_df, on='Term', how='inner')
    df.sort_values(by='count', inplace=True, ascending=False)

    # 创建自定义颜色映射
    cmap = cm.get_cmap(color_bar)  # 使用红绿渐变色阶
    norm = Colors.Normalize(vmin=df[column].min(),
                            vmax=df[column].max())
    colors = [cmap(norm(value)) for value in df[column][0:n_count]]

    # 获取所有的go_type
    go_types = df['go_type'].unique()

    for i, go_type in enumerate(go_types):
        # 绘制柱状图
        fig, ax = plt.subplots()
        # 获取该go_type的数据
        sub_df = df[df['go_type'] == go_type]
        colors = cmap(norm(sub_df[column]))
        bars = ax.barh(sub_df['Term'][0:n_count], sub_df['count'][0:n_count], color=colors)
        # 创建颜色条
        sm = cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])  # 必须设置一个空数组

        cbar = plt.colorbar(sm)
        cbar.set_label(column)

        # 设置图形属性
        ax.set_xlabel('Count')
        ax.set_ylabel('Term')
        ax.set_title(f"{pre}_{anno_type}_barplot_{go_type}")

        ofname = os.path.join(enrichment_out, format_filename(f"{pre}_{anno_type}_barplot_{go_type}") + ".pdf")
        plt.savefig(ofname, bbox_inches='tight', pad_inches=0.2)
        plt.clf()

    # 计算气泡图的大小和位置
    for i, go_type in enumerate(go_types):
        # 获取该go_type的数据
        sub_df = df[df['go_type'] == go_type]
        x = (sub_df["Odds Ratio"][0:n_count])
        y = sub_df['Term'][0:n_count]
        sizes = (sub_df["count"][0:n_count] * 20)
        colors = [cmap(norm(value)) for value in sub_df[column][0:n_count]]

        # 绘制气泡图
        fig, ax = plt.subplots()
        scatter = ax.scatter(x, y, s=sizes, c=colors)

        # 创建颜色条
        sm = cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])  # 必须设置一个空数组

        cbar = plt.colorbar(sm)
        cbar.set_label(column)

        # 设置图形属性
        ax.set_xlabel('Odds Ratio')
        ax.set_title(f"{pre}_{anno_type}_dotplot_{go_type}")

        # 增加表示气泡大小的图例
        handles, labels = scatter.legend_elements(prop="sizes", alpha=0.1, func=lambda x: x / 20)
        legend2 = ax.legend(handles, labels, loc="upper right", title="Count", bbox_to_anchor=(1.4, 1), borderaxespad=0)

        ofname = os.path.join(enrichment_out, format_filename(f"{pre}_{anno_type}_dotplot_{go_type}") + ".pdf")
        plt.savefig(ofname, bbox_inches='tight', pad_inches=0.2)
        plt.clf()
    # ofname=os.path.join(enrichment_out, format_filename(f"{pre}_{anno_type}_dotplot") + ".pdf")
    # fig = ax.figure  # 获取Figure对象
    # fig.savefig(ofname, bbox_inches='tight', pad_inches=0)


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
    # config = "./config.py"

    # 获得config文件的内容
    global pipeline
    pipeline = False
    if config != "":
        # 获得config文件的内容
        import sys
        global cfg
        if not os.path.exists(config):
            raise Exception(f"配置文件{config}不存在，请检查路径")
        try:
            module_name = os.path.basename(config).split('.')[0]
            cfg = __import__(module_name)
        except:
            sys.path.append(os.path.dirname(config))
            cfg = __import__(module_name)

        global config_prefix
        config_prefix = os.path.splitext(config)[0]
        faa_args = cfg.faa_args
        n_count = faa_args.get('n_count')
        color_bar = faa_args.get('color_bar')
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
    merged_df = pd.merge(deg_df, blast_df, on="query_id", how="inner")

    annotation_df = pd.read_csv(annotation_file, sep='\t')
    annotation_df = annotation_df.rename(columns={annotation_df.columns[0]: 'entry'})
    annotation_df = pd.merge(blast_df, annotation_df, on="entry", how="inner")
    global go_type_df
    go_type_df = get_go_type(annotation_df)

    # 获得差异表达基因的各类注释信息
    deg_out = os.path.join(output, deg_out)
    deg_out_df = get_out(annotation_df, merged_df, deg_out)

    # 获得所有可鉴别转录本的各类注释信息
    all_out = os.path.join(output, all_out)
    all_out_df = get_out(annotation_df, blast_df, all_out)

    # 获得背景信息：blast后的到的
    global background
    background = blast_df["query_id"].tolist()
    background = list(set(background))

    # 富集分析与绘图
    enrichment_out = os.path.join(output, enrichment_out)
    if not os.path.exists(enrichment_out):
        os.makedirs(enrichment_out)

    # 上调基因：
    up_gene_list = merged_df[merged_df["logFC"] > 0]["query_id"].tolist()
    up_gene_list = list(set(up_gene_list))
    # 下调基因：
    down_gene_list = merged_df[merged_df["logFC"] < 0]["query_id"].tolist()
    down_gene_list = list(set(down_gene_list))

    # pre = "up"
    # gene_list = up_gene_list

    # 先进行GO富集分析
    anno_type = "Gene Ontology (GO)"
    gene_sets = get_gene2anno(annotation_df, anno_type, outdir=output)
    if gene_sets != {}:
        logger.info(f"开始进行{anno_type}富集分析\n")
        # 上调：
        draw_write_go("up", up_gene_list, gene_sets, anno_type, n_count, color_bar, enrichment_out)
        # 下调
        draw_write_go("down", down_gene_list, gene_sets, anno_type, n_count, color_bar, enrichment_out)
        logger.info(f"{anno_type}富集分析完成！\n")

    # 进行kegg通路富集分析
    anno_type = "Pathway"
    gene_sets = get_gene2anno(annotation_df, anno_type, outdir=output)
    if gene_sets != {}:
        logger.info(f"开始进行{anno_type}富集分析\n")
        annotation_df[anno_type] = annotation_df[anno_type].str.replace('PATHWAY: ', '')
        # 上调：
        draw_write("up", up_gene_list, gene_sets, anno_type, n_count, color_bar, enrichment_out)
        # 下调
        draw_write("down", down_gene_list, gene_sets, anno_type, n_count, color_bar, enrichment_out)
        logger.info(f"{anno_type}富集分析完成！\n")

    # 其他的分类
    if len(annotation_df.columns) > 14:
        for i in range(14, len(annotation_df.columns)):
            anno_type = annotation_df.columns[i]
            gene_sets = get_gene2anno(annotation_df, anno_type, outdir=output)
            if gene_sets != {}:
                logger.info(f"开始进行{anno_type}富集分析\n")
                # 上调：
                draw_write("up", up_gene_list, gene_sets, anno_type, n_count, color_bar, enrichment_out)
                # 下调
                draw_write("down", down_gene_list, gene_sets, anno_type, n_count, color_bar, enrichment_out)
                logger.info(f"{anno_type}富集分析完成！\n")

    # 读取原始 pkl 文件中的列表数据
    if pipeline:
        config_pkl = f'{config_prefix}.pkl'
        with open(config_pkl, "rb") as f:
            data = pickle.load(f)

        # 在原始列表中添加新的字符串元素
        enrichment_out_file = os.listdir(enrichment_out)  # 获取out下的所有文件和文件夹
        res_dict = {"all_out": all_out, "deg_out": deg_out, "enrichment_out": enrichment_out_file}
        data["annotation_analysis_files"] = res_dict

        # 将原始列表和新的字符串一起写入新的 pkl 文件中
        with open(config_pkl, "wb") as f:
            pickle.dump(data, f)
    logger.info(f"特征注释分析完成！\n")


if __name__ == '__main__':
    main()
