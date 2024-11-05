#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@File    :   Transcript_identification.py    
@Contact :   2454888366@qq.com

@Modify Time      @Author    @Version    @Desciption
------------      -------    --------    -----------
2023/3/25 上午10:50   skychou      1.0         None
'''
import logging as lg
import os
import pickle
import re
import shutil
import click


# 判断序列数据库文件是否存在
def db_exist(blast_db):
    dblist = [blast_db + ".pin", blast_db + ".nsq", blast_db + ".nhr"]
    flag = True
    for i in dblist:
        if not os.path.exists(i):
            return False
    return flag

    # 简要统计信息
    # 读取blast结果文件


def get_res_file(blast_res_file, stat, summary, transcript_out):
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
    all_identified_transcripts = ["query_id\tentry"]
    # 统计结果文件
    blast_stat_list = ["#统计结果如下：\nquery_id\thits_found\n"]
    blast_stat_dict = dict()
    for i in every_blast:
        query_id = i[0].split()[2]
        if len(i) <= 3:
            hits = i[2].split()[1]
        else:
            hits = i[3].split()[1]
            tempquery_id = i[4].split()[0]
            tempentry = i[4].split()[1]
            all_identified_transcripts.append(f"{tempquery_id}\t{tempentry}")
        if not blast_stat_dict.__contains__(hits):
            blast_stat_dict[hits] = [query_id]
        else:
            blast_stat_dict[hits].append(query_id)
        blast_stat_list.append(f"{query_id}\t{hits}\n")
        sort_dict_list = sorted(blast_stat_dict.items(), key=lambda x: x[0])

    blast_summary_list = ["摘要信息如下：\n"]
    for i in sort_dict_list:
        blast_summary_list.append(f"\n{i[0]} hits found的共有{str(len(i[1]))}条\n")
        blast_summary_list.append(f"{','.join(i[1])}\n")

    with open(stat, "w+") as f:
        f.writelines(blast_stat_list)
    with open(summary, "w+") as f:
        f.writelines(blast_summary_list)
    with open(transcript_out, "w") as f:
        f.write("\n".join(all_identified_transcripts))


@click.command()
@click.option('--input_file', '-i', type=click.Path(), help='组装转录本序列文件', prompt=True)
@click.option('--protein_file', '-p', type=click.Path(), help='已知蛋白序列文件', default="")
@click.option('--summary_file', '-s', type=click.Path(), help='summary表格文件', default="")
@click.option("--query", "-q",
              help='检索关键词，例如"(taxonomy_id:4751) AND (reviewed:true)"\n详细请查看手册：https://www.uniprot.org/help/query-fields',
              default="")
@click.option("--addition_fields", "-af", default="",
              help='检索该物种信息的额外条目(列)，用逗号分割，默认为:"accession", "xref_geneid","gene_names", "organism_id", "organism_name", "cc_pathway", "go_id", "go"')
@click.option('--output', '-o', type=str, help='结果输出文件(绝对路径)', default="./transcript_indentification_output/")
@click.option('--outfmt', '-of', type=str, default="7", help='结果输出格式(默认为7)')
@click.option('--evalue', '-e', default=10 ** -5, help='evalue值(默认为10e-5)')
@click.option('--cpu', type=str, default="1", help='CPU核数(默认为1)')
@click.option("--config", "-c", type=click.Path(), default="",
              help="配置文件目录")
def blast_main_op(input_file, protein_file, summary_file, query, addition_fields, output, outfmt, evalue, cpu, config):
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
        logger.error(f"输出路径{output_dir}不存在！\n")
        raise Exception(f"输出路径{output_dir}不存在！")

    if not output_dir.endswith("/"):
        output_dir = output_dir + "/"
    output = output_dir

    # 没有给定handler日志级别，将用logger的级别
    fileHandler = lg.FileHandler(filename=f"{output_dir}transcript_identification.log", mode="w", encoding="UTF-8")
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

    # 参数检验
    if not os.path.exists(input_file):
        logger.error(f"输入文件{input_file}不存在！\n")
        raise Exception(f"输入文件{input_file}不存在！")

    if protein_file and summary_file:
        # 如果是protein_file和summary_file文件存在，则无需执行子程序Uniport_query.py
        # 复制文件到目标目录
        p_basename = os.path.basename(protein_file)
        shutil.copy(protein_file, os.path.join(output, p_basename))
        protein_file = os.path.join(output, p_basename)
        fasta_output = protein_file

        s_basename = os.path.basename(summary_file)
        shutil.copy(summary_file, os.path.join(output, s_basename))
        summary_output = os.path.join(output, s_basename)


        if pipeline:
            # 读取原始 pkl 文件中的列表数据
            config_pkl = f'{config_prefix}.pkl'
            with open(config_pkl, "rb") as f:
                data = pickle.load(f)

            # 在原始列表中添加新的字符串元素
            tempdict = {"fasta_file": os.path.join(fasta_output),
                        "summary_file": os.path.join(summary_output)}
            data["uniport_query_files"] = tempdict

            # 将原始列表和新的字符串一起写入新的 pkl 文件中
            with open(config_pkl, "wb") as f:
                pickle.dump(data, f)

    elif query:
        # 如果是query参数，则执行子程序Uniport_query.py
        if addition_fields!="":
            cmd = f"python3 Uniprot_query.py --query '{query}' --addition_fields '{addition_fields}' --output_dir {output} --config {config}"
        else:
            cmd = f"python3 Uniprot_query.py --query '{query}' --output_dir {output} --config {config}"
        rc = os.system(cmd)
        # 判断是否执行成功
        if rc == 0:
            logger.info(f"Uniprot_query.py执行成功！\n")
        else:
            logger.error(f"Uniprot_query.py执行失败！\n")
            raise Exception(f"Uniprot_query.py执行失败！")
        if pipeline:
            config_pkl = f'{config_prefix}.pkl'
            if os.path.exists(config_pkl):
                with open(config_pkl, 'rb') as f:
                    files_dict = pickle.load(f)
                    protein_file = files_dict["uniport_query_files"]["fasta_file"]
            else:
                logger.error(f"没有找到pickle文件{config_pkl}！\n")
                raise Exception(f"没有找到pickle文件{config_pkl}！")
        else:
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

            file_prefix = format_filename(query)
            protein_file = file_prefix + ".fasta.gz"
            protein_file = os.path.join(output, protein_file)
    else:
        logger.error(f"请输入必要参数已知蛋白序列文件(--protein_file)或检索关键词(--query)!\n")
        raise Exception("请输入必要参数已知蛋白序列文件(--protein_file)或检索关键词(--query)!")

    if os.path.exists(protein_file):
        if protein_file.endswith("gz"):
            fasta_file = protein_file[:-3]
            cmd = f"gzip -c -d {protein_file} > {fasta_file}"
            rc = os.system(cmd)
            # 判断是否执行成功
            if rc == 0:
                logger.info(f"gzip执行成功！\n")
            else:
                logger.error(f"gzip执行失败！\n")
                raise Exception(f"gzip执行失败！")
            # print(cmd)
        else:
            fasta_file = protein_file

    # 执行makeblastdb
    blast_db_dir = os.path.join(output, "blastdb")
    if not os.path.exists(blast_db_dir):
        os.makedirs(blast_db_dir)
    blastdb = "RSNAP_db"
    cmd = f"makeblastdb -in {fasta_file} -dbtype prot -title {blastdb} -out {os.path.join(blast_db_dir, blastdb)} -parse_seqids"
    rc = os.system(cmd)
    # 判断是否执行成功
    if rc == 0:
        logger.info(f"makeblastdb执行成功！\n")
    else:
        logger.error(f"makeblastdb执行失败！\n")
        raise Exception(f"makeblastdb执行失败！")
    filename = os.path.basename(fasta_file)  # file.txt
    prefix = os.path.splitext(filename)[0]  # file
    blast_result = os.path.join(output, f"{prefix}.outfmt7")
    if not os.path.exists(blast_result):
        cmd = f"blastx -query {input_file} -db {os.path.join(blast_db_dir, blastdb)} -out {blast_result} -outfmt {outfmt} -evalue {evalue} -num_threads {cpu} -max_target_seqs 1"
        rc = os.system(cmd)
        # 判断是否执行成功
        if rc == 0:
            logger.info(f"blastx执行成功！\n")
        else:
            logger.error(f"blastx执行失败！\n")
            raise Exception(f"blastx执行失败！")
    else:
        logger.info(f"blast_result文件:{blast_result}已存在！\n")

    stat_file = os.path.join(output, "blastx.stat")
    summary_file = os.path.join(output, "blastx.summary")
    transcript_out = os.path.join(output, "all_identified_transcripts.tsv")
    get_res_file(blast_result, stat_file, summary_file, transcript_out)
    logger.info(f"Transcript_identification.py执行完成！\n")

    # 读取原始 pkl 文件中的列表数据
    config_pkl = f'{config_prefix}.pkl'
    with open(config_pkl, "rb") as f:
        data = pickle.load(f)

    # 在原始列表中添加新的字符串元素
    tempdict = {"blast_result": blast_result,
                "stat_file": stat_file,
                "summary_file": summary_file,
                "transcript_out": transcript_out,
                "fasta_file": fasta_file}
    data["transcript_indentification_files"] = tempdict

    # 将原始列表和新的字符串一起写入新的 pkl 文件中
    with open(config_pkl, "wb") as f:
        pickle.dump(data, f)


if __name__ == '__main__':
    blast_main_op()
