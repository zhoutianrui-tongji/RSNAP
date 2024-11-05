#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@File    :   RSNAP.py    
@Contact :   2454888366@qq.com

@Modify Time      @Author    @Version    @Desciption
------------      -------    --------    -----------
2023/4/11 上午11:30   skychou      1.0         None
'''

# if not os.path.exists("./.first/"):
# 检查是否是第一次执行
#     import os
#     import sys
#     import subprocess
#
#     if sys.version_info < (3, 8):
#         raise Exception("python({})版本过低,请更新至3.7以上版本！".format(sys.version))
#     current_dir = os.path.dirname(os.path.abspath(__file__))
#     os.chdir(current_dir)
#     # 尝试导入依赖
#     try:
#         import matplotlib
#
#         matplotlib.use('Agg')
#
#         import gseapy
#         import numpy
#         import glob
#         import zipfile
#         import shutil
#         import logging as lg
#         import os
#         import shutil
#         import pickle
#         import click
#         import rpy2
#         import pandas
#         import re
#         import time
#         import urllib
#         import requests
#     except ImportError as e:
#         # 获取缺少的包名
#         package = e.name
#         # 调用pip命令来安装
#         print(f"缺少python包:{package}，正在安装！如果安装失败请手动安装：pip install {package}")
#         subprocess.run(["pip", "install", package])
#     # 检查R以及R包
#     try:
#         from rpy2.robjects import r
#
#         # 运行R的命令
#         output = subprocess.check_output(['R', '--version'])
#         version_string = output.decode('utf-8').split()[2]  # 从输出中获取版本字符串
#
#         # 提取主版本号和次版本号
#         major_version, minor_version = map(int, version_string.split('.')[:2])
#     except:
#         raise Exception("rpy2导入失败，未能连接到R！")
#     if major_version < 3 or (major_version == 3 and minor_version < 8):
#         raise Exception("R({})版本过低,请更新至3.8以上版本！".format(major + "." + minor))
#     from rpy2.robjects.packages import importr
#
#     utils = importr('utils')
#
#
#     def extract_versions(package_data):
#         return dict(zip(
#             package_data.rx(True, 'Package'),  # get Package column
#             package_data.rx(True, 'Version')  # get Version column
#         ))
#
#
#     installed_packages = extract_versions(utils.installed_packages())
#     r_packages = {'limma': '3.50.3', 'edgeR': '3.36.0',
#                   'pheatmap': '1.0.12', 'data.table': '1.14.8'}  # 'DESeq2': '1.34.0',
#     for packnames in r_packages.items():
#         if not (packnames[0] in installed_packages):
#             print(f"缺少R包:{packnames[0]}，正在安装！如果安装失败请手动安装：install.packages('{packnames[0]}')")
#             utils.chooseCRANmirror(ind=1)  # select a CRAN mirror
#             utils.install_packages(packnames[0])  # install the package
#         else:
#             need_version = ".".join(packnames[1].split(".")[:-1])
#             packages_version = ".".join(installed_packages[packnames[0]].split(".")[:-1])
#             if float(packages_version) < float(need_version):
#                 print(
#                     f"R包:{packnames[0]}版本过低，正在升级！如果升级失败请手动升级：library(remotes)install_version('{packnames[0]}', '{packnames[1]}')")
#                 utils.chooseCRANmirror(ind=1)  # select a CRAN mirror
#                 utils.install_packages(packnames[0])  # install the package
#
#     # 检查和安装linux命令
#     linux_packages = ["fastqc", "bowtie2", "pigz", "gzip"]
#     for linux_package in linux_packages:
#         code = subprocess.run(["which", linux_package]).returncode  # 检查fastqc命令是否存在
#         if code != 0:
#             print(f"缺少必备程序:{linux_package}，正在安装！如果安装失败请手动安装：apt-get install {linux_package}")
#             raise Exception(f"缺少必备程序:{linux_package}!")
#             # subprocess.run(["apt", "install", "fastqc"])  # 安装fastqc命令
#             # subprocess.run(["apt", "update"])  # 更新软件源信息
#             # subprocess.run(["apt", "upgrade"])  # 更新所有安装的ubuntu命令
#     some_packages = ["TrimmomaticPE", "trimmomatic"]
#     for linux_package in some_packages:
#         code = subprocess.run(["which", linux_package]).returncode  # 检查fastqc命令是否存在
#         if code == 0:
#             break  #
#     else:
#         print(f"缺少必备程序:{linux_package}，正在安装！如果安装失败请手动安装：apt-get install {linux_package}")
#         raise Exception(f"缺少必备程序:{linux_package}!")
#     extend_software = ["Trinity", "RSEM"]
#     for ext in extend_software:
#         if ext == "Trinity":
#             linux_package = "Trinity"
#             code = subprocess.run(["which", linux_package]).returncode  # 检查fastqc命令是否存在
#             if code != 0:
#                 if subprocess.run(["which", "conda"]).returncode:
#                     subprocess.run(["conda", "install", "-c", "bioconda", "trinity"])
#                     print(
#                         f"缺少必备程序:{ext}，正在使用conda安装！如果安装失败请手动安装：conda install -c bioconda trinity")
#                 else:
#                     print(f"缺少必备程序:{ext}，conda不存在！正在下载编译安装！如果安装失败请手动安装\n："
#                           f"wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/Trinity-v2.15.1/trinityrnaseq-v2.15.1.FULL.tar.gz\n"
#                           f"tar -xvzf trinityrnaseq-v2.15.1.FULL.tar.gz\n"
#                           f"cd trinityrnaseq-v2.15.1\n"
#                           f"make\n"
#                           f"make install")
#                     cmds = [
#                         f"wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/Trinity-v2.15.1/trinityrnaseq-v2.15.1.FULL.tar.gz",
#                         f"tar -xvzf trinityrnaseq-v2.15.1.FULL.tar.gz",
#                         f"cd trinityrnaseq-v2.15.1",
#                         f"make",
#                         f"make install"
#                     ]
#                     for install_cmd in cmds:
#                         try:
#                             os.system(install_cmd)
#                         except Exception as e:
#                             raise Exception(f"安装{linux_package}失败！请手动安装！")
#         if ext == "RSEM":
#             linux_package = "rsem-prepare-reference"
#             code = subprocess.run(["which", linux_package]).returncode  # 检查fastqc命令是否存在
#             if code != 0:
#                 raise Exception(f"安装{ext}失败！请手动安装！")
#     # 创建一个隐藏的文件夹
#     else:
#         os.makedirs("./.first/")

# 导入必须的包
import os
import click
import pickle
import logging as lg
import sys
import subprocess
import re


@click.command()
@click.option("--config", "-c", type=click.Path(), default="./config.py",
              help="配置文件目录，默认为./config.py")
@click.option("--check", is_flag=True,default=False,
              help="是否调用cheke.py子程序进行依赖检查，默认为否")
def main(config,check):
    import os
    if check:
        print("开始检查依赖...")
        rc = os.system(os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "./check.py")))
        # 判断是否执行成功
        if rc == 0:
            print(f"依赖检查成功！\n")
        else:
            print(f"依赖检查失败！\n")
            raise Exception(f"依赖检查失败！")



    # config = "./config.py"  ##测试
    current_dir = os.path.dirname(os.path.abspath(__file__))
    # 获得config文件的内容
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

    # 获得用户输入
    user_args = cfg.user_args
    samples = user_args.get("samples")
    groups = user_args.get("groups")
    layout = user_args.get("layout").lower()
    format = user_args.get("format").lower()
    output_dir = user_args.get("output_dir")
    #转换为绝对路径
    output_dir = os.path.abspath(output_dir)
    query = user_args.get("query")
    addition_fields = user_args.get("addition_fields")
    compare = user_args.get("compare")

    # out 目录检查并处理
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

    # 设置日志文件
    # 记录器
    global logger
    logger = lg.getLogger("mylog")
    logger.setLevel(lg.DEBUG)

    # 处理器
    consleHandler = lg.StreamHandler()
    consleHandler.setLevel(lg.DEBUG)

    # 没有给定handler日志级别，将用logger的级别
    global log_file
    log_file = f"{output_dir}RSNAP.log"
    fileHandler = lg.FileHandler(filename=log_file, mode="w+", encoding="UTF-8")
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

    logger.info(f"开始执行RSNAP主程序......\n")

    os.chdir(current_dir)
    ################################################################
    # 执行Quality_Control.py

    # 获得Quality_Control参数
    qc_args = cfg.qc_args
    not_cut = qc_args.get("not_cut")
    threshold = qc_args.get("threshold")
    threads = qc_args.get("threads")
    software = qc_args.get("software")
    all_samples = []
    for sample in samples.values():
        all_samples.append(f"{sample[0]},{sample[1]}")
    input = "%".join(all_samples)
    output = os.path.join(output_dir, "qc_fastqc_output/")
    qc_fastqc_output_dir = output
    script = qc_args.get("script")

    flag = False
    # 检查结果文件是否存在：
    if os.path.exists(output):
        files = os.listdir(output)  # 获取目录下的所有文件和文件夹

        def get_pre(filename):
            filename = os.path.basename(filename)  # file.txt
            prefix = os.path.splitext(filename)[0]  # file
            suffix = os.path.splitext(filename)[1]
            if suffix == ".gz":
                prefix = os.path.splitext(prefix)[0]
            return prefix

        prefix_list = [get_pre(sample) for sample in all_samples]
        # 是否所有的文件都存在

        for prefix in prefix_list:
            if prefix not in "#".join(files):
                break
        else:
            flag = True

    pkl_flag = os.path.exists(f'{config_prefix}.pkl')

    if (not flag) or (not pkl_flag):
        if not_cut:
            not_cut = " --not_cut "
        else:
            not_cut = " "
        # 拼接命令行
        qc_cmd = f"./Quality_Control.py " \
                 f"-i '{input}' -o '{output}' --format {format}{not_cut}--threads {threads} " \
                 f"--threshold {threshold} --software {software} --config {config}"
        logger.info(f"开始执行Quality_Control.py......\n")
        logger.info(f"命令如下：{qc_cmd}")
        try:
            # 使用subprocess.run来执行命令，并捕获输出
            result = subprocess.run(qc_cmd, shell=True, check=True, text=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            # 如果命令执行失败，记录错误输出并抛出异常
            logger.error(f"Quality_Control.py执行失败！错误输出：{e.stderr}")
            raise Exception(f"Quality_Control.py执行失败！错误输出：{e.stderr}")
        else:
            # 如果命令执行成功，记录标准输出
            logger.info("Quality_Control.py执行成功！输出结果：\n")
    else:
        logger.info(f"Quality_Control.py结果文件已经存在：{output}\n将执行后续步骤！\n")
    files = os.listdir(output)  # 获取目录下的所有文件和文件夹
    qc_files = [os.path.join(output, file) for file in files]
    ################################################################
    # 执行Assemble_trinity.py
    # 获得Assemble_trinity参数
    at_args = cfg.at_args
    result_file = at_args.get("result_file")
    summary_file = at_args.get("summary_file")
    max_memory = at_args.get("max_memory")
    cpu = at_args.get("cpu")
    output = os.path.join(output_dir, "assemble_trinity_output/")
    assemble_trinity_output_dir = output
    script = at_args.get("script")

    flag = False
    if os.path.exists(output):
        files = os.listdir(output)  # 获取目录下的所有文件和文件夹
        # 检查文件是否存在
        if (result_file in "#".join(files)) and (summary_file in "#".join(files)):
            flag = True

    # 读取trim_files
    if os.path.exists(f'{config_prefix}.pkl'):
        with open(f'{config_prefix}.pkl', 'rb') as f:
            data = pickle.load(f)
            trim_files = data['trim_files']

    if not flag:
        if script == "./Assemble_trinity.py":
            if format == "fasta":
                tempformat = "fa"
            elif format == "fastq":
                tempformat = "fq"
            else:
                logger.error(f"format参数设置错误！{format}/n")
                raise Exception(f"format参数设置错误！{format}")

            # 如果是双端测序
            if layout == "paired":
                left = []
                right = []
                for sample in trim_files.values():
                    left.append(sample[0])
                    right.append(sample[1])
                left = ",".join(left)
                right = ",".join(right)
                at_cmd = f"./Assemble_trinity.py" \
                         f" --format {tempformat} --paired -l {left} -r {right} --output {output} --summary_file {summary_file}" \
                         f" --max_memory {max_memory} --cpu {cpu} --config {config}"
            elif layout == "single":
                input = []
                for sample in trim_files.values():
                    input.append(sample[0])
                input = ",".join(input)
                at_cmd = f"./Assemble_trinity.py" \
                         f" --format {tempformat} --input {input} --output {output} --summary_file {summary_file}" \
                         f" --max_memory {max_memory} --cpu {cpu} --config {config}"
            else:
                logger.error(f"layout参数设置错误！{layout}/n")
                raise Exception(f"layout参数设置错误！{layout}")
            logger.info(f"开始执行Assemble_trinity.py......\n")
            logger.info(f"命令如下：{at_cmd}")
            try:
                # 使用subprocess.run来执行命令，并捕获输出
                result = subprocess.run(at_cmd, shell=True, check=True, text=True, stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE)
            except subprocess.CalledProcessError as e:
                # 如果命令执行失败，记录错误输出并抛出异常
                logger.error(f"Assemble_trinity.py执行失败！错误输出：{e.stderr}")
                raise Exception(f"Assemble_trinity.py执行失败！错误输出：{e.stderr}")
            else:
                # 如果命令执行成功，记录标准输出
                logger.info("Assemble_trinity.py执行成功\n")
        else:
            ###后续拓展在这里写
            pass
    else:
        logger.info(f"Assemble_trinity.py结果文件已经存在：{output}\n将执行后续步骤！\n")
    files = os.listdir(output)  # 获取目录下的所有文件和文件夹
    at_files = [os.path.join(output, file) for file in files]

    # 读取assemble_files
    if script == "./Assemble_trinity.py":
        with open(f'{config_prefix}.pkl', 'rb') as f:
            data = pickle.load(f)
            assemble_file = data['assemble_files']["assembly_file"]
            summary_file = data['assemble_files']["summary_file"]
    else:
        pass
        # 自己的读取方法，就是要获得组装结果文件assemble_file和总结文件summary_file的绝对路径

    ################################################################
    # 执行Expression_evaluation_rsem.py
    # 获得Expression_evaluation_rsem.py参数
    eer_args = cfg.eer_args
    threads = eer_args.get("threads")
    output = os.path.join(output_dir, "expression_rsem_output/")
    expression_rsem_output_dir = output
    script = eer_args.get("script")

    flag = False
    if os.path.exists(output):
        files = os.listdir(output)  # 获取目录下的所有文件和文件夹
        # 检查文件是否存在
        if ("genes.expected_count.matrix" in "#".join(files)):
            flag = True

    if not flag:
        if script == "./Expression_evaluation_rsem.py":
            input = ",".join(all_samples)
            # 如果是双端测序
            if layout == "paired":
                eer_cmd = f"./Expression_evaluation_rsem.py --paireds {input} --reference {assemble_file} --output {output} --config {config} --threads {threads}"
            else:
                eer_cmd = f"./Expression_evaluation_rsem.py --singles {input} --reference {assemble_file} --output {output} --config {config} --threads {threads}"
            logger.info(f"开始执行Expression_evaluation_rsem.py......\n")
            logger.info(f"命令如下：{eer_cmd}")
            try:
                # 使用subprocess.run来执行命令，并捕获输出
                result = subprocess.run(eer_cmd, shell=True, check=True, text=True, stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE)
            except subprocess.CalledProcessError as e:
                # 如果命令执行失败，记录错误输出并抛出异常
                logger.error(f"Expression_evaluation_rsem.py执行失败！错误输出：{e.stderr}")
                raise Exception(f"Expression_evaluation_rsem.py执行失败！错误输出：{e.stderr}")
            else:
                # 如果命令执行成功，记录标准输出
                logger.info("Expression_evaluation_rsem.py执行成功\n")
        else:
            pass
        ###后续拓展在这里写
    else:
        logger.info(f"Expression_evaluation_rsem.py结果已经存在：{output}\n将执行后续步骤！\n")
    files = os.listdir(output)  # 获取目录下的所有文件和文件夹
    eer_files = [os.path.join(output, file) for file in files]

    if script == "./Expression_evaluation_rsem.py":
        # 读取expression_matrix
        with open(f'{config_prefix}.pkl', 'rb') as f:
            data = pickle.load(f)
            expression_matrix = data["expression_matrix"]
    else:
        pass
        # 获得表达矩阵expression_matrix的绝对路径

    ################################################################
    # Transcript_identification.py
    # 获得Transcript_identification.py参数
    ti_args = cfg.ti_args
    cpu = ti_args.get("cpu")
    evalue = ti_args.get("evalue")
    output = os.path.join(output_dir, "transcript_indentification_output/")
    transcript_indentification_output_dir = output
    script = ti_args.get("script")
    protein_file = user_args.get("protein_file")
    uq_summary_file = user_args.get("uq_summary_file")


    if addition_fields!="":
        ti_cmd = f"./Transcript_identification.py" \
                 f" --input_file {assemble_file} --protein_file {protein_file} --summary_file {uq_summary_file} --query '{query}' --addition_fields {addition_fields} --output {output}" \
                 f" --evalue {evalue} --cpu {cpu} --config {config}"
    else:
        ti_cmd = f"./Transcript_identification.py" \
                 f" --input_file {assemble_file} --protein_file {protein_file} --summary_file {uq_summary_file} --query '{query}' --output {output}" \
                 f" --evalue {evalue} --cpu {cpu} --config {config}"

    logger.info(f"开始执行Transcript_identification.py......\n")
    logger.info(f"命令如下：{ti_cmd}")
    try:
        # 使用subprocess.run来执行命令，并捕获输出
        result = subprocess.run(ti_cmd, shell=True, check=True, text=True, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        # 如果命令执行失败，记录错误输出并抛出异常
        logger.error(f"Transcript_identification.py执行失败！错误输出：{e.stderr}")
        raise Exception(f"Transcript_identification.py执行失败！错误输出：{e.stderr}")
    else:
        # 如果命令执行成功，记录标准输出
        logger.info("Transcript_identification.py执行成功\n")
    files = os.listdir(output)  # 获取目录下的所有文件和文件夹
    ti_files = [os.path.join(output, file) for file in files]

    for compared in compare:
        ################################################################
        # Differential_expression_gene_identification.py
        # 获得Differential_expression_gene_identification.py参数
        degi_args = cfg.degi_args
        software = degi_args.get("software")
        foldchange = degi_args.get("foldchange")
        padj = degi_args.get("padj")
        name_pre = "_vs_".join(compared)
        output = os.path.join(output_dir, f"{name_pre}_DEG_identification_output/")
        DEG_identification_output_dir = output
        script = degi_args.get("script")

        flag = False
        if os.path.exists(output):
            files = os.listdir(output)  # 获取目录下的所有文件和文件夹
            # 检查文件是否存在
            if ("sigDEG.tsv" in "#".join(files)) and (summary_file in "#".join(files)):
                flag = True

        if not flag:
            if script == "./Differential_expression_gene_identification.py":
                degi_cmd = f'''./Differential_expression_gene_identification.py --input_file {expression_matrix} --output_dir {output} --software {software} --foldchange {foldchange} --padj {padj} --compare "{compared}" --config {config}'''

                logger.info(f"开始执行Differential_expression_gene_identification.py({name_pre})......\n")
                logger.info(f"命令如下：{degi_cmd}")
                try:
                    # 使用subprocess.run来执行命令，并捕获输出
                    result = subprocess.run(degi_cmd, shell=True, check=True, text=True, stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE)
                except subprocess.CalledProcessError as e:
                    # 如果命令执行失败，记录错误输出并抛出异常
                    logger.error(f"Differential_expression_gene_identification.py执行失败！错误输出：{e.stderr}")
                    raise Exception(f"Differential_expression_gene_identification.py执行失败！错误输出：{e.stderr}")
                else:
                    # 如果命令执行成功，记录标准输出
                    logger.info("Differential_expression_gene_identification.py执行成功\n")
            else:
                pass
                ###后续拓展在这里写
        else:
            logger.info(f"Differential_expression_gene_identification.py结果文件({name_pre})已经存在：{output}\n将执行后续步骤！\n")
        files = os.listdir(output)  # 获取目录下的所有文件和文件夹
        degi_files = [os.path.join(output, file) for file in files]

        if script == "./Differential_expression_gene_identification.py":
            # 读取deg_identification_files
            with open(f'{config_prefix}.pkl', 'rb') as f:
                data = pickle.load(f)
                deg_identification_files = data["deg_identification_files"]
        else:
            pass
            # 获得差异分析文件deg_identification_files的绝对路径字典
            # 例如：
            # # 获取out下的所有文件和文件夹
            # all_items = os.listdir(output)
            # # 过滤掉文件夹，只保留所有的文件
            # deg_identification_files = [os.path.join(output,i) for i in all_items]

        ################################################################
        # Feature_annotation_analysis.py
        # Feature_annotation_analysis.py参数
        faa_args = cfg.faa_args
        all_out = faa_args.get("all_out")
        deg_out = faa_args.get("deg_out")
        enrichment_out = faa_args.get("enrichment_out")
        output = os.path.join(output_dir, f"{name_pre}_feature_annotation_analysis_output/")
        Feature_annotation_analysis_dir = output
        script = faa_args.get("script")

        # 读取uniport_query_files和transcript_indentification_files
        with open(f'{config_prefix}.pkl', 'rb') as f:
            data = pickle.load(f)
            uniport_query_files = data["uniport_query_files"]
            annotation_file = uniport_query_files["summary_file"]
            transcript_indentification_files = data["transcript_indentification_files"]
            fasta_file = transcript_indentification_files["fasta_file"]
            blast_list = transcript_indentification_files["transcript_out"]
            blastx_summary_file = transcript_indentification_files["summary_file"]
            deg_list = [i for i in deg_identification_files if "sigDEG.tsv" in i][0]

        flag = False
        if os.path.exists(output):
            files = os.listdir(output)  # 获取目录下的所有文件和文件夹
            # 检查文件是否存在
            if (enrichment_out in "#".join(files)) and (all_out in "#".join(files)) and (deg_out in "#".join(files)):
                flag = True

        if not flag:
            faa_cmd = f"./Feature_annotation_analysis.py" \
                      f" --annotation_file {annotation_file} --blast_list {blast_list} --deg_list {deg_list}" \
                      f" --output {output} --all_out {all_out} --deg_out {deg_out} --enrichment_out {enrichment_out}" \
                      f" --config {config}"
            logger.info(f"开始执行Feature_annotation_analysis.py({name_pre})......\n")
            logger.info(f"命令如下：{faa_cmd}")
            try:
                # 使用subprocess.run来执行命令，并捕获输出
                result = subprocess.run(faa_cmd, shell=True, check=True, text=True, stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE)
            except subprocess.CalledProcessError as e:
                # 如果命令执行失败，记录错误输出并抛出异常
                logger.error(f"Feature_annotation_analysis.py执行失败！错误输出：{e.stderr}")
                raise Exception(f"Feature_annotation_analysis.py执行失败！错误输出：{e.stderr}")
            else:
                # 如果命令执行成功，记录标准输出
                logger.info("Feature_annotation_analysis.py执行成功\n")
        else:
            logger.info(f"Feature_annotation_analysis.py结果文件已经存在({name_pre})：{output}\n")

        with open(f'{config_prefix}.pkl', 'rb') as f:
            data = pickle.load(f)
            annotation_analysis_files = data['annotation_analysis_files']
            deg_annotation_file = annotation_analysis_files["deg_out"]
            transcripts_annotation_file = annotation_analysis_files["all_out"]
        files = os.listdir(output)  # 获取目录下的所有文件和文件夹
        faa_files = [os.path.join(output, file) for file in files]
        faa_photo_files = os.listdir(os.path.join(output, enrichment_out))
        faa_photo_files = [os.path.join(os.path.join(output, enrichment_out), file) for file in faa_photo_files]
        ###################################################
        # 生成html报告
        from jinja2 import Environment, FileSystemLoader
        import os
        import pandas as pd
        # 获得总表
        # 表达矩阵
        eer_matrix_df = pd.read_csv(os.path.join(expression_rsem_output_dir, "genes.expected_count.matrix"), sep='\t')
        # 使用正则表达式从列名中提取 S1, S2, S3 等名称
        new_cols = eer_matrix_df.columns.str.extract(r'output/(.*)\.genes')
        # 将列名修改为新的名称
        eer_matrix_df.columns = new_cols.squeeze().values
        eer_matrix_df = eer_matrix_df.rename(columns={eer_matrix_df.columns[0]: 'ID'})

        # 差异分析结果
        degi_allGenes_out = [i for i in degi_files if "allGenes" in i][0]
        degi_allGenes_out_df = pd.read_csv(degi_allGenes_out, sep=' ')
        degi_allGenes_out_df = degi_allGenes_out_df.reset_index()
        degi_allGenes_out_df = degi_allGenes_out_df.rename(columns={'index': 'ID'})

        res_df = pd.merge(eer_matrix_df, degi_allGenes_out_df, on='ID', how='outer')

        # 读取blast结果outfmt7文件
        blast_result = [i for i in ti_files if "outfmt" in i][0]
        # 读取blast结果outfmt7文件
        blast_result_df = pd.read_csv(blast_result, comment='#', sep='\s+', header=None)
        # 将数据框的各列进行重命名
        blast_result_df.columns = ['ID', 'subject_id', 'perc_identity', 'aln_length', 'mismatch', 'gap_open', 'q_start',
                                   'q_end', 's_start', 's_end', 'evalue', 'bit_score']
        res_df = pd.merge(res_df, blast_result_df, on='ID', how='outer')

        all_anno = [i for i in faa_files if "all_anno" in i][0]
        all_anno_df = pd.read_csv(all_anno, sep='\t')
        all_anno_df = all_anno_df.rename(columns={all_anno_df.columns[0]: 'ID'})

        res_df = pd.merge(res_df, all_anno_df, on='ID', how='outer')
        res_df.to_csv(os.path.join(output_dir, f'all_data({name_pre}).csv'), index=False)

        logger.info(f'({name_pre})总表保存成功！')

        # assemble_file = "./assemble_trinity_output/assemblies.fa"  ##测试
        # summary_file = "./assemble_trinity_output/summary.txt"  ##测试

        # Define the template environment
        env = Environment(loader=FileSystemLoader(current_dir))
        template = env.get_template('report_template.html')

        # Define the data variables
        # ##测试
        # qc_fastqc_output_dir = 'qc_fastqc_output'
        # assemble_trinity_output_dir = 'assemble_trinity_output'
        # expression_rsem_output_dir = 'expression_rsem_output'
        # DEG_identification_output_dir = 'DEG_identification_output'
        # transcript_indentification_output_dir = 'transcript_indentification_output'
        # Feature_annotation_analysis_dir = 'feature_annotation_analysis_output'

        trinity_assembly = "".join(open(assemble_file).readlines()[:100])
        trinity_summary = "".join(open(summary_file).readlines()[:100])

        rsem_matrix = "".join(
            open(os.path.join(expression_rsem_output_dir, 'genes.expected_count.matrix')).readlines()[:100])

        for file in degi_files:
            if re.findall(r"FDR.*logFC", file):
                deg_table_file = file
                break
        deg_table = "".join(
            open(deg_table_file).readlines()[:100])
        # deg_table = "".join(
        #     open(os.path.join(DEG_identification_output_dir, 'DESeq2_sigDEGs_FDR0.05_logFC1.out')).readlines()[:100])  ##测试

        uniport_summary = "".join(open(annotation_file).readlines()[:100])
        # uniport_summary = "".join(open(
        #     os.path.join(transcript_indentification_output_dir, 'taxonomy_id_4751_and_reviewed_true.summary')).readlines()[
        #                           :100]) #测试

        uniport_sequences = "".join(open(fasta_file).readlines()[:100])
        # uniport_sequences = "".join(open(
        #     os.path.join(transcript_indentification_output_dir, 'taxonomy_id_4751_and_reviewed_true.fasta')).readlines()[
        #                             :100]) #测试

        identified_transcripts = "".join(open(blast_list).readlines()[:100])
        # identified_transcripts = "".join(open(
        #     os.path.join(transcript_indentification_output_dir, 'all_identified_transcripts.tsv')).readlines()[:100])
        blastx_summary = "".join(open(blastx_summary_file).readlines()[:100])
        # blastx_summary = "".join(
        #     open(os.path.join(transcript_indentification_output_dir, 'blastx.summary')).readlines()[:100])

        deg_annotation = "".join(open(deg_annotation_file).readlines()[:100])
        # deg_annotation = "".join(open(os.path.join(Feature_annotation_analysis_dir, 'deg_anno.tsv')).readlines()[:100])

        transcripts_annotation = "".join(open(transcripts_annotation_file).readlines()[:100])
        # transcripts_annotation = "".join(
        #     open(os.path.join(Feature_annotation_analysis_dir, 'all_anno.tsv')).readlines()[:100])

        qc_describe = f'在这一部分中，使用FastQC对样本进行质量控制，然后根据质控结果，使用{qc_args.get("software")}对原始测序文件进行裁剪，用于后续分析。'
        qc_reference = f'''<span class="ref">Simon Andrews. FastQC: A quality control tool for high throughput sequence data. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/. Accessed on 27 May 2023. Version 0.12.0.</span>
                        <br>
                        <span class="ref">Bolger, A. M., Lohse, M., & Usadel, B. (2014). Trimmomatic: A flexible trimmer for Illumina Sequence Data. Bioinformatics, 30(15), 2114-2120. doi:10.1093/bioinformatics/btu170</span>
    '''
        templist = []
        for file in qc_files:
            if file.endswith(".html") or "summary" in file:
                templist.append(file)
        templist.sort()
        tempres = ""
        for i in range(1, len(templist), 2):
            summary = templist[i]
            html_file = templist[i - 1]
            with open(summary, "r", encoding='utf-8') as f:
                lines = f.readlines()
            sample_name = lines[0].strip().split(":")[1][:-11]
            Q20 = lines[1].strip().split(":")[1].strip()
            Q30 = lines[2].strip().split(":")[1].strip()
            PASS = lines[3].strip().split("!")[1].strip()
            header = lines[4].strip().split(":")[1]
            tail = lines[5].strip().split(":")[1]
            tempres += f'''
              <tr>
        <td>{sample_name}</td>
        <td>{Q20}</td>
        <td>{Q30}</td>
        <td>{PASS}</td>
        <td>{header}</td>
        <td>{tail}</td>
        <td><a href="{html_file.replace(output_dir, './')}" target=“_blank”>{os.path.split(html_file)[1]}</a></td>
      </tr>'''
            fastqc_res = f'''<table>
              <tr>
                <th>样本名</th>
                <th>Q20</th>
                <th>Q30</th>
                <th>通过检查</th>
                <th>开头剪切碱基数</th>
                <th>尾部剪切碱基数</th>
                <th>详细报告</th>
              </tr>
              {tempres}
            </table>'''
        # for file in qc_files:
        #     if file.endswith(".html"):
        #         fastqc_res = fastqc_res + f'''
        #                 <tr>
        #                 <td><a href="{file.replace(output_dir, './')}" target=“_blank”>{os.path.split(file)[1]}</a></td>
        #                 </tr>'''
        # fastqc_summary = ""
        # for file in qc_files:
        #     if "summary" in file:
        #         fastqc_summary += f"""  <p><strong>{os.path.split(file)[1]}</strong></p>
        #                                 <pre>{"".join(open(file).readlines())}</pre>"""



        logfile = [i for i in at_files if i.endswith("log")][0]
    #     assemblies_res = f'''
    #                     <tr>
    #                     <td><a href="{assemble_file.replace(output_dir, './')}" target=“_blank”>组装结果文件：{os.path.split(assemble_file)[1]}</a></td>
    #                     </tr>
    #                     <tr>
    #                     <td><a href="{summary_file.replace(output_dir, './')}" target=“_blank”>组装摘要文件{os.path.split(summary_file)[1]}</a></td>
    #                     </tr>
    # '''
        with open(summary_file, 'r') as f:
            lines = f.readlines()
        tolal_gene = lines[5].split()[-1]
        tolal_transcript = lines[6].split()[-1]
        at_describe = f'''在这一部分，使用Trinity对进行裁剪后的测序数据进行从头组装，并且调用TrinityStats.pl生成组装摘要信息。得到了<a href="{assemble_file.replace(output_dir, './')}" target=“_blank”>组装结果文件：{os.path.split(assemble_file)[1]}</a>
        和<a href="{summary_file.replace(output_dir, './')}" target=“_blank”>组装摘要文件{os.path.split(summary_file)[1]}</a>
        组装后的到的基因有{tolal_gene}个，转录本共有{tolal_transcript}个。'''
        at_reference = f'''<span class="ref">Grabherr, M. G., Haas, B. J., Yassour, M., Levin, J. Z., Thompson, D. A., Amit, I., … & Regev, A. (2011). Full-length transcriptome assembly from RNA-Seq data without a reference genome. Nature biotechnology, 29(7), 644-652. doi:10.1038/nbt.1883</span>
        '''


        expected_matrix = os.path.join(expression_rsem_output_dir, 'genes.expected_count.matrix')
        eer_describe = f'''在这一部分，使用RSEM软件对RNA-seq数据进行表达量定量分析，并生成<a href="{expected_matrix.replace(output_dir, './')}" target=“_blank”>表达矩阵文件：{os.path.split(expected_matrix)[1]}</a>'''
        eer_reference = f'''<span class="ref">Li, B., & Dewey, C. N. (2011). RSEM: accurate transcript quantification from RNA-Seq data with or without a reference genome. BMC bioinformatics, 12(1), 323. doi:10.1186/1471-2105-12-323</span>
        '''

        degi_describe = f'''在这一部分，使用python:rpy2包使python与R交互，用{degi_args.get("software")}包进行差异分析，并使用ggplot2包和pheatmap包绘制热图和火山图等图进行可视化。得到的文件有'''

        def get_line(file):
            with open(file, 'r') as f:
                lines = f.readlines()
                return len(lines)-1

        for file in degi_files:
            if file.endswith("out"):
                if "down" in file:
                    degi_describe = degi_describe + f'''<a href="{file.replace(output_dir, './')}" target=“_blank”>所有下调基因差异分析结果：{os.path.split(file)[1]}</a>,共有{get_line(file)}个。'''
                elif "up" in file:
                    degi_describe = degi_describe + f'''<a href="{file.replace(output_dir, './')}" target=“_blank”>所有上调基因差异分析结果：{os.path.split(file)[1]}</a>,共有{get_line(file)}个。'''
                elif "allGenes" in file:
                    degi_describe = degi_describe + f'''<a href="{file.replace(output_dir, './')}" target=“_blank”>所有基因差异分析结果：{os.path.split(file)[1]}</a>,共有{get_line(file)}个。'''
                elif "sigDEGs" in file and "FDR" in file and "logFC" in file:
                    degi_describe = degi_describe + f'''<a href="{file.replace(output_dir, './')}" target=“_blank”>按FDR(阈值：{padj})和logFC(阈值：{foldchange})筛选的差异基因差异分析结果：{os.path.split(file)[1]}</a>,共有{get_line(file)}个。'''
                else:
                    degi_describe = degi_describe + f'''<a href="{file.replace(output_dir, './')}" target=“_blank”>按FDR(阈值：{padj})筛选的差异基因差异分析结果：{os.path.split(file)[1]}</a>,共有{get_line(file)}个。'''

        degi_photo = ""
        # index = 0
        # for file in degi_files:
        #     if file.endswith("pdf"):
        #         if index == 0:
        #             degi_photo += f"""
        #     <div class="carousel-item active pdf-container">
        #      <embed class="d-block" src="{file.replace(output_dir, './')}#toolbar=0" type='application/pdf'>
        #     </div>
        #     """
        #         else:
        #             degi_photo += f"""
        #             <div class="carousel-item pdf-container">
        #               <embed class="d-block" src="{file.replace(output_dir, './')}#toolbar=0" type='application/pdf'>
        #             </div>
        #             """
        #         index += 1
        degi_files.sort()
        count = 0
        for file in degi_files:
            if file.endswith("pdf"):
                count += 1
                if "volcano" in file:
                    degi_photo += f'''
                    <br>
                        <div class="pdf-container">
                           <embed class="d-block" src="{file.replace(output_dir, './')}#toolbar=0" type='application/pdf'>
                        </div>
                        <p class="centered-text">图{count}：火山图</p> '''
                if "heatmap" in file:
                    degi_photo += f'''
                                    <br>
                        <div class="pdf-container">
                           <embed class="d-block" src="{file.replace(output_dir, './')}#toolbar=0" type='application/pdf'>
                        </div>
                        <p class="centered-text">图{count}：热图</p> '''
                if "scatter" in file:
                    degi_photo += f'''
                    <br>
                        <div class="pdf-container">
                           <embed class="d-block" src="{file.replace(output_dir, './')}#toolbar=0" type='application/pdf'>
                        </div>
                        <p class="centered-text">图{count}：scatter map 散点图</p> '''
                if "MDS" in file:
                    degi_photo += f'''
                    <br>
                        <div class="pdf-container">
                           <embed class="d-block" src="{file.replace(output_dir, './')}#toolbar=0" type='application/pdf'>
                        </div>
                        <p class="centered-text">图{count}：多维缩放（MDS）聚类图</p> '''
                if "clustering_tree" in file:
                    degi_photo += f'''
                    <br>
                        <div class="pdf-container">
                           <embed class="d-block" src="{file.replace(output_dir, './')}#toolbar=0" type='application/pdf'>
                        </div>
                        <p class="centered-text">图{count}：层次聚类树状图</p> '''
        degi_reference = f'''<span class="ref">Li, B. (2013). rpy2: A Python interface to R. https://rpy2.github.io/. Accessed on 27 May 2023. Version 3.5.13.</span>
                            <br><span class="ref">Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol. 2014;15(12):550. doi:10.1186/s13059-014-0550-8</span>
                            <br><span class="ref">Robinson MD, McCarthy DJ and Smyth GK (2010). edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics 26, 139-140</span>
                            <br><span class="ref">Wickham H. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016</span>
                            <br><span class="ref">Kolde R (2019). pheatmap: Pretty Heatmaps. R package version 1.0.12. https://CRAN.R-project.org/package=pheatmap2</span>
                            <br><span class="ref">Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W, Smyth GK (2015). limma powers differential expression analyses for RNA-sequencing and microarray studies. Nucleic Acids Research 43(7):e472</span>
        '''

        ti_describe = f"在这一部分，使用blastx工具将组装转录本序列文件与已知蛋白序列文件或检索关键词对应的UniProt数据库中的蛋白序列进行比对，然后根据比对结果识别转录本并进行统计和摘要分析。得到的结果文件有："
        for file in ti_files:
            if file.endswith("fasta"):
                ti_describe = ti_describe + f'''<a href="{file.replace(output_dir, './')}" target=“_blank”>Uniport下载的测序文件：{os.path.split(file)[1]};</a>'''
            if file.endswith("fasta.gz"):
                ti_describe = ti_describe + f'''<a href="{file.replace(output_dir, './')}" target=“_blank”>Uniport下载的测序文件(压缩)：{os.path.split(file)[1]};</a>'''
            if file.endswith("outfmt7"):
                ti_describe = ti_describe + f'''<a href="{file.replace(output_dir, './')}" target=“_blank”>Blastx结果文件：{os.path.split(file)[1]};</a>'''
            if file.endswith("summary"):
                if "blastx.summary" in file:
                    ti_describe = ti_describe + f'''<a href="{file.replace(output_dir, './')}" target=“_blank”>Blastx结果总结：{os.path.split(file)[1]};</a>'''
                else:
                    ti_describe = ti_describe + f'''<a href="{file.replace(output_dir, './')}" target=“_blank”>Uniport下载的注释信息文件：{os.path.split(file)[1]};</a>'''
            if file.endswith("stat"):
                ti_describe = ti_describe + f'''<a href="{file.replace(output_dir, './')}" target=“_blank”>Blastx结果统计信息：{os.path.split(file)[1]};</a>'''
                df = pd.read_table(file, sep='\t', header=1)
                counts = df['hits_found'].value_counts()
                num_zero = counts[0]
                num_nonzero = len(df) - num_zero
            if file.endswith("tsv"):
                ti_describe = ti_describe + f'''<a href="{file.replace(output_dir, './')}" target=“_blank”>所有鉴别到的转录本：{os.path.split(file)[1]};</a>'''
                all_ti_num = get_line(file)

        ti_describe += f"在Blast结果中，0 hits found 的转录本有{num_zero}条,非 0 hits found的 转录本有{num_nonzero}条。所有鉴别到的转录本共有{all_ti_num}条。"
        ti_reference = f'''<span class="ref">Altschul, S.F., Gish, W., Miller, W., Myers, E.W., Lipman, D.J. (1990) “Basic local alignment search tool.” J. Mol. Biol. 215:403-410. PubMed</span>
        <br><span class="ref">The UniProt Consortium. UniProt: a worldwide hub of protein knowledge Nucleic Acids Res. 47:D506-515 (2019)</span>
        '''



        faa_describe = f"在这一部分，使用python:gseapy包进行富集分析，并且进行可视化。得到的结果文件有："
        for file in faa_files:
            if "all_anno" in file:
                faa_describe = faa_describe + f'''<a href="{file.replace(output_dir, './')}" target=“_blank”>所有鉴定到的转录本的注释信息：{os.path.split(file)[1]};</a>'''
            if "deg_anno" in file:
                faa_describe = faa_describe + f'''<a href="{file.replace(output_dir, './')}" target=“_blank”>差异基因的注释信息：{os.path.split(file)[1]};</a>'''
            if "GO" in file:
                faa_describe = faa_describe + f'''<a href="{file.replace(output_dir, './')}" target=“_blank”>GO富集分析结果表：{os.path.split(file)[1]};</a>'''
            if "Pathway" in file:
                faa_describe = faa_describe + f'''<a href="{file.replace(output_dir, './')}" target=“_blank”>KEGG通路富集分析结果表：{os.path.split(file)[1]};</a>'''

        faa_photo = ""
        # index = 0
        # for file in faa_photo_files:
        #     if file.endswith("pdf"):
        #         if index == 0:
        #             faa_photo += f"""
        #         <div class="carousel-item active pdf-container">
        #          <embed class="d-block" src="{file.replace(output_dir, './')}#toolbar=0" type='application/pdf'>
        #         </div>
        #         """
        #         else:
        #             faa_photo += f"""
        #                 <div class="carousel-item pdf-container">
        #                   <embed class="d-block" src="{file.replace(output_dir, './')}#toolbar=0" type='application/pdf'>
        #                 </div>
        #                 """
        #         index += 1
        count = 0
        faa_photo_files.sort()
        for file in faa_photo_files:
            if file.endswith("pdf"):
                count += 1
                if "down_gene_ontology_go_barplot" in file:
                    for temp in ["cc", "bp", "mf"]:
                        if temp in file:
                            temp = temp.upper()
                            break
                    faa_photo += f'''
                        <br>
                        <div class="pdf-container">
                           <embed class="d-block" src="{file.replace(output_dir, './')}#toolbar=0" type='application/pdf'>
                        </div>
                        <p class="centered-text">图{count}:下调基因GO({temp})富集分析柱状图</p> '''
                if "down_gene_ontology_go_dotplot" in file:
                    for temp in ["cc", "bp", "mf"]:
                        if temp in file:
                            temp = temp.upper()
                            break
                    faa_photo += f'''
                    <br>
                        <div class="pdf-container">
                           <embed class="d-block" src="{file.replace(output_dir, './')}#toolbar=0" type='application/pdf'>
                        </div>
                        <p class="centered-text">图{count}:下调基因GO({temp})富集分析气泡图</p> '''
                if "down_pathway_barplot" in file:
                    faa_photo += f'''
                    <br>
                        <div class="pdf-container">
                           <embed class="d-block" src="{file.replace(output_dir, './')}#toolbar=0" type='application/pdf'>
                        </div>
                        <p class="centered-text">图{count}:下调基因KEGG通路富集分析柱状图</p> '''
                if "down_pathway_dotplot" in file:
                    faa_photo += f'''
                    <br>
                        <div class="pdf-container">
                           <embed class="d-block" src="{file.replace(output_dir, './')}#toolbar=0" type='application/pdf'>
                        </div>
                        <p class="centered-text">图{count}:下调基因KEGG通路富集分析气泡图</p> '''
                if "up_gene_ontology_go_barplot" in file:
                    for temp in ["cc", "bp", "mf"]:
                        if temp in file:
                            temp = temp.upper()
                            break
                    faa_photo += f'''
                    <br>
                        <div class="pdf-container">
                           <embed class="d-block" src="{file.replace(output_dir, './')}#toolbar=0" type='application/pdf'>
                        </div>
                        <p class="centered-text">图{count}:上调基因GO({temp})富集分析柱状图</p> '''
                if "up_gene_ontology_go_dotplot" in file:
                    for temp in ["cc", "bp", "mf"]:
                        if temp in file:
                            temp = temp.upper()
                            break
                    faa_photo += f'''
                    <br>
                        <div class="pdf-container">
                           <embed class="d-block" src="{file.replace(output_dir, './')}#toolbar=0" type='application/pdf'>
                        </div>
                        <p class="centered-text">图{count}:上调基因GO({temp})富集分析气泡图</p> '''
                if "up_pathway_barplot" in file:
                    faa_photo += f'''
                    <br>
                        <div class="pdf-container">
                           <embed class="d-block" src="{file.replace(output_dir, './')}#toolbar=0" type='application/pdf'>
                        </div>
                        <p class="centered-text">图{count}:上调基因KEGG通路富集分析柱状图</p> '''
                if "up_pathway_dotplot" in file:
                    faa_photo += f'''
                    <br>
                        <div class="pdf-container">
                           <embed class="d-block" src="{file.replace(output_dir, './')}#toolbar=0" type='application/pdf'>
                        </div>
                        <p class="centered-text">图{count}: 上调基因KEGG通路富集分析气泡图</p> '''

        faa_reference = f'''<span class="ref">Zhuoqing Fang, Xinyuan Liu, Gary Peltz, GSEApy: a comprehensive package for performing gene set enrichment analysis in Python, Bioinformatics, 2022;, btac757, https://doi.org/10.1093/bioinformatics/btac757</span>
                            <br><span class="ref">Kanehisa, M. and Goto, S.; KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic Acids Res. 28, 27-30 (2000). PubMed</span>
                            <br><span class="ref">Ashburner, M., Ball, C.A., Blake, J.A., Botstein, D., Butler, H., Cherry, J.M., Davis, A.P., Dolinski, K., Dwight, S.S., Eppig, J.T., Harris, M.A., Hill, D.P., Issel-Tarver, L., Kasarskis, A., Lewis, S., Matese, J.C., Richardson, J.E., Ringwald, M., Rubin, G.M. and Sherlock, G.; Gene ontology: tool for the unification of biology. The Gene Ontology Consortium. Nat Genet. 25:25-29 (2000). PubMed</span>
    '''
        # Render the template
        html = template.render(
            qc_describe=qc_describe,
            qc_reference=qc_reference,
            at_describe=at_describe,
            at_reference=at_reference,
            eer_describe=eer_describe,
            eer_reference=eer_reference,
            degi_describe=degi_describe,
            degi_reference=degi_reference,
            ti_describe=ti_describe,
            ti_reference=ti_reference,
            faa_describe=faa_describe,
            faa_reference=faa_reference,
            fastqc_res=fastqc_res,
            trinity_assembly=trinity_assembly,
            trinity_summary=trinity_summary,
            rsem_matrix=rsem_matrix,
            degi_photo=degi_photo,
            # deg_table=deg_table,
            # uniport_summary=uniport_summary,
            # uniport_sequences=uniport_sequences,
            identified_transcripts=identified_transcripts,
            blastx_summary=blastx_summary,
            # deg_annotation=deg_annotation,
            # transcripts_annotation=transcripts_annotation,
            faa_photo=faa_photo
        )

        # Write the HTML report to a file
        with open(os.path.join(output_dir, f'report({name_pre}).html'), 'w') as f:
            f.write(html)

    logger.info("RSNAP执行成功！")


if __name__ == '__main__':
    main()
