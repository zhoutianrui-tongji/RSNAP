#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@File    :   Assemble_trinity.py
@Contact :   2454888366@qq.com

@Modify Time      @Author    @Version    @Desciption
------------      -------    --------    -----------
2/26/2023 9:48 AM   skychou      1.0         None
'''
import logging as lg
import os
import pickle
import shutil
import subprocess

# 导入click包
import click


# 定义一个函数，用于执行trinity rsem组装
def run_trinity(format, input_files, left, right, output_dir, paired, result_file, max_memory, cpu):
    # 根据单端或双端数据，设置不同的参数
    if not paired:
        cmd = f"Trinity --seqType {format} --single {input_files} --max_memory {max_memory}G --CPU {cpu} --output {output_dir[:-1]}"
    else:
        cmd = f"Trinity --seqType {format} --left {left} --right {right} --max_memory {max_memory}G --CPU {cpu} --output {output_dir[:-1]}"
    # 执行命令
    if not os.path.exists(os.path.join(output_dir, result_file)):
        logger.info(f"开始执行Trinity......\n")
        rc = os.system(cmd)
        # 判断是否执行成功
        if rc == 0:
            logger.info(f"Trinity执行成功！\n")
        else:
            logger.error(f"Trinity执行失败！\n")
            raise Exception(f"Trinity执行失败！")

        # 返回组装结果文件的路径

        # 移动结果文件并且改名
        # 还需测试
        '''python ./Assemble_trinity.py --paired -l ./test/Con-1_1.fq.gz -r ./test/Con-1_2.fq.gz -o ./assemble_trinity_output/Con-1'''
        fasta_file = f"{output_dir[:-1]}.Trinity.fasta"
        if os.path.exists(fasta_file):
            gene_trans_map = fasta_file + ".gene_trans_map"
            # shutil.move(fasta_file, output_dir)
            # shutil.move(gene_trans_map, output_dir)
            os.rename(gene_trans_map,
                      os.path.join(output_dir, "Trinity.fasta.gene_trans_map"))
        else:
            fasta_file = f"{output_dir}Trinity.fasta"
        os.rename(fasta_file, os.path.join(output_dir, result_file))
    else:
        logger.info(f"Trinity结果文件：{result_file}已存在！将直接进行后续步骤！\n")
    return os.path.join(output_dir, result_file)


# 定义一个函数，用于提取组装摘要信息，并输出到指定文件中
def extract_assembly_summary(assembly_file, summary_file):
    if subprocess.getoutput("which TrinityStats.pl") == "":
        # 先找到TrinityStats.pl
        # result = subprocess.getstatusoutput("find / -name TrinityStats.pl | grep -m 1 'TrinityStats.pl'")
        result = subprocess.getstatusoutput("locate TrinityStats.pl")
        lines = result[1].split()
        pl_script = i
        for i in lines:
            if "TrinityStats.pl" in i:
                break
        else:
            logger.error(
                f"没有找到TrinityStats.pl！请确认是否安装trinityrnaseq！\n")
            raise Exception(
                f"没有找到TrinityStats.pl！请确认是否安装trinityrnaseq！")
    else:
        pl_script = "TrinityStats.pl"
    # 使用trinity自带的脚本，生成组装摘要信息
    cmd = f"{pl_script} {assembly_file} > {summary_file}"
    # 执行命令
    rc = os.system(cmd)
    # 判断是否执行成功
    if rc == 0:
        logger.info(f"提取组装摘要信息成功！\n")
    else:
        logger.error(f"提取组装摘要信息失败！\n")
        raise Exception(f"提取组装摘要信息失败！")

    # 也可以自己编写
    # 提取组装摘要信息，并按照指定格式输出到指定文件中


# 定义一个click命令行接口，用于接收用户输入的参数
@click.command()
@click.option("--format", "-fm", type=click.Choice(["fq", "fa"]), default="fa",
              help="测序数据文件格式,默认为'fasta'")
@click.option("--input", "-i", type=click.Path(), help="单端测序数据存放目录或文件（多个文件使用逗号隔开）")
@click.option("--left", "-l", type=click.Path(), help="双端测序中(左段)数据存放目录或文件（多个文件使用逗号隔开）")
@click.option("--right", "-r", type=click.Path(), help="双端测序中(右段)数据存放目录或文件（多个文件使用逗号隔开）")
@click.option("--output", "-o", type=click.Path(), default="./assemble_trinity_output", help="组装结果文件")
@click.option("--paired", "-p", is_flag=True, help="是否为双端测序数据，默认为否")
@click.option("--result_file", "-f", type=click.Path(), default="assemblies.fa", help="组装结果文件名")
@click.option("--summary_file", "-s", type=click.Path(), default="summary.txt", help="组装摘要结果文件")
@click.option("--max_memory", "-m", default="6", help="最大内存占用数（G）（默认为6G）")
@click.option("--cpu", "-c", default="4", help="cup核数，默认为4")
@click.option("--config", "-cfg", type=click.Path(), default="",
              help="配置文件目录")
def main(format, input, left, right, output, paired, result_file, summary_file, max_memory, cpu, config):
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
        raise Exception("out dir does not exist!输出路径不存在！")

    if not output_dir.endswith("/"):
        output_dir = output_dir + "/"
    output = output_dir

    # 没有给定handler日志级别，将用logger的级别
    fileHandler = lg.FileHandler(filename=f"{output_dir}assemble_trinity.log", mode="w", encoding="UTF-8")
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

    # 检查输入参数是否合法
    if (not input and not paired) or (not left or not right):
        logger.error(f"请输入必要的参数！\n")
        raise Exception(f"请输入必要的参数！")

    # 之前以为是一个一个组装
    # 如果输入多个文件，使用逗号分割，并转换为列表
    # input_files = input
    # if not paired:
    #     if "," in input:
    #         input_files = input.split(",")
    #     else:
    #         input_files = input
    # else:
    #     if "," in left and "," in right:
    #         left = left.split(",")
    #         right = right.split(",")

    # 调用run_trinity函数，进行组装，并获取组装结果文件路径

    assembly_file = run_trinity(format, input, left, right, output, paired, result_file, max_memory, cpu)
    # 调用extract_assembly_summary函数，提取组装摘要信息，并输出到指定文件中
    summary_file = os.path.join(output, summary_file)
    extract_assembly_summary(assembly_file, summary_file)

    if pipeline:
        # 读取原始 pkl 文件中的列表数据
        config_pkl = f'{config_prefix}.pkl'
        with open(config_pkl, "rb") as f:
            data = pickle.load(f)

        # 在原始列表中添加新的字符串元素
        data["assemble_files"] = {"assembly_file": assembly_file, "summary_file": summary_file}

        # 将原始列表和新的字符串一起写入新的 pkl 文件中
        with open(config_pkl, "wb") as f:
            pickle.dump(data, f)


# 如果是直接运行该脚本，则执行main函数
if __name__ == "__main__":
    main()
