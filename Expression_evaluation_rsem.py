#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@File    :   Expression_evaluation_rsem.py
@Contact :   2454888366@qq.com

@Modify Time      @Author    @Version    @Desciption
------------      -------    --------    -----------
3/01/2023 3:28 PM   skychou      1.0         None
'''

# 导入logging包
import logging as lg
# 导入os包
import os

# 导入click包
import click

import pickle

# 获得最大前缀公共子字符串
def longest_common_prefix(str_list):
    # 检查列表是否为空
    if not str_list:
        return ""
    # 对列表中的所有字符串进行排序
    str_list.sort()
    # 初始化前缀为空字符串
    prefix = ""
    # 遍历第一个和最后一个字符串的每个字符
    for i in range(min(len(str_list[0]), len(str_list[-1]))):
        # 如果字符相同，将它们添加到前缀中
        if str_list[0][i] == str_list[-1][i]:
            prefix += str_list[0][i]
        # 否则，跳出循环
        else:
            break
    # 返回前缀
    return prefix

# 定义一个命令行函数，使用@click装饰器
@click.command()
# 添加一些选项参数，如单端/双端、测序数据文件、参考转录本文件等
@click.option('--singles', '-s', help='单端测序数据文件列表（用,分割）', default="")
@click.option('--paireds', '-p', help='双端测序数据文件列表（,分割）（两端文件依次排列）', default="")
@click.option('--reference', '-r', required=True, help='参考转录本文件')
@click.option('--output', '-o', type=click.Path(), default='./expression_rsem_output', help='输出目录')
@click.option("--config", "-c", type=click.Path(), default="./config.py",
              help="配置文件目录，默认为./config.py")
@click.option("--threads", "-t", default="8", help="线程数，默认为8")
# 定义函数主体，使用click.echo打印信息
def rsem_quantify(singles, paireds, reference, output, config, threads):
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

    # 从文件中读取数据并反序列化
    # 获得分组和样本信息
    global trim_files
    if os.path.exists(f'{config_prefix}.pkl'):
        with open(f'{config_prefix}.pkl', 'rb') as f:
            data = pickle.load(f)
            trim_files = data['trim_files']
            is_pipeline = True
    else:
        is_pipeline = False

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
    global log_file
    log_file = f"{output_dir}expression_evaluation_rsem.log"
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
    if is_pipeline:
        if cfg.user_args["layout"] == "paired":
            templist = []
            for i in trim_files.values():
                templist.extend(i)
            paireds = ",".join(templist)
        if cfg.user_args["layout"] == "single":
            templist = []
            for i in trim_files.values():
                templist.extend(i)
            singles = ",".join(templist)

    # 检查输入参数是否合法
    if singles and paireds:
        logger.error(f"错误：不能同时指定单端和双端测序数据文件！\n")
        raise Exception(f"错误：不能同时指定单端和双端测序数据文件！")
    if not singles and not paireds:
        logger.error(f"错误：必须指定单端或双端测序数据文件！\n")
        raise Exception(f"错误：必须指定单端或双端测序数据文件！")
    if not os.path.exists(reference):
        logger.error(f"错误：参考文件不存在！\n")
        raise Exception(f"错误：参考文件不存在！")

    # 打印开始信息
    if os.path.exists(f"{output}genes.expected_count.matrix"):
        logger.info(f"{output}genes.expected_count.matrix文件已经存在，将执行后续步骤！\n")
        return
    else:
        logger.info(f"开始运行trinity rsem...\n")

    # 构造rsem-prepare-reference命令，使用reference作为输入，output作为输出目录
    '''
    rsem-prepare-reference ./test/Con-1.fasta ./sample/sample1_output/Trinity --bowtie2
    '''
    '''
    rsem-prepare-reference ./test/Con-1.fasta ./mytest/Con-1/Trinity --bowtie2 -p 8
    '''
    prepare_cmd = f'rsem-prepare-reference {reference} {output}Trinity --bowtie2'
    # 调用系统命令执行rsem-prepare-reference
    rc = os.system(prepare_cmd)
    # 判断是否执行成功
    if rc == 0:
        logger.info(f"rsem-prepare-reference命令执行成功！\n")
    else:
        logger.error(f"rsem-prepare-reference命令执行失败！\n")
        raise Exception(f"rsem-prepare-reference命令执行失败！")
    # print(prepare_cmd)
    # 之前的错误想法，以为要分开处理
    # # 拆解reference参数
    # try:
    #     references = references.split(",")
    # except Exception:
    #     click.echo('错误：references必须使用;分割！')
    #     return

    # prefix_list = []
    # for reference in references:
    #     # 获取前缀名
    #     filename = os.path.basename(reference)  # file.txt
    #     prefix = os.path.splitext(filename)[0]  # file
    #     prefix_list.append(prefix)
    #
    #     if not os.path.exists(f"{output}{prefix}"):
    #         os.makedirs(f"{output}{prefix}")
    #     else:
    #         print("路径已经存在，请确认是否覆盖！")
    #
    #     prepare_cmd = f'rsem-prepare-reference {reference} {output}{prefix}/Trinity --bowtie2 -p {threads}'
    #     # 调用系统命令执行rsem-prepare-reference
    #     os.system(prepare_cmd)
    #     # print(prepare_cmd)

    # 根据单端或双端构造rsem-calculate-expression命令，使用single或paired作为输入，output作为输出目录，correction作为修正参数
    '''rsem-calculate-expression --paired-end ./test/Con-1_1.fq.gz ./test/Con-1_2.fq.gz --bowtie2 -p 8 --no-bam-output ./sample/sample1_output/Trinity ./sample/sample1_output/Con-1'''
    '''rsem-calculate-expression --paired-end ./test/Con-1_1.fq.gz ./test/Con-1_2.fq.gz --bowtie2 -p 8 --no-bam-output ./mytest/Con-1/Trinity ./mytest/Con-1/Con-1'''
    #记录结果文件
    reslist = []
    if singles:
        # 拆解singles参数
        try:
            singles = singles.split(",")
        except Exception:
            logger.error(f"错误：singles必须使用,分割！\n")
            raise Exception(f"错误：singles必须使用,分割！")

        # 之前的错误想法，以为要分开处理
        # for prefix, single in enumerate(singles):
        #     calculate_cmd = f'rsem-calculate-expression --single-cell-prior {single} --bowtie2 -p {threads} --no-bam-output {output}{prefix}/Trinity {output}{prefix}/{prefix}'
        #     # print(calculate_cmd)
        #     os.system(calculate_cmd)

        for single in singles:
            filename = os.path.basename(single)  # file.txt
            if not is_pipeline:
                prefix = os.path.splitext(filename)[0]  # file
            else:
                inverse_dict = {tuple(v): k for k, v in trim_files.items()}
                key = (single,)
                prefix = inverse_dict[key]
            calculate_cmd = f'rsem-calculate-expression --single-cell-prior {single} --bowtie2 -p {threads} --no-bam-output {output}Trinity {output}{prefix}'
            # print(calculate_cmd)
            reslist.append(f"{output}{prefix}.genes.results")
            rc = os.system(calculate_cmd)
            # 判断是否执行成功
            if rc == 0:
                logger.info(f"rsem-calculate-expression命令执行成功！\n")
            else:
                logger.error(f"rsem-calculate-expression命令执行失败！\n")
                raise Exception(f"rsem-calculate-expression命令执行失败！")

    else:
        # 拆解paireds参数
        try:
            paireds = paireds.split(",")
        except Exception:
            logger.error(f"错误：paireds必须使用,分割！\n")
            raise Exception(f"错误：paireds必须使用,分割！")

        paireds = [(paireds[i], paireds[i + 1]) for i in range(0, len(paireds), 2)]

        # 之前的错误想法，以为要分开处理
        # for n, paired in enumerate(paireds):
        #     prefix = prefix_list[n]
        #     calculate_cmd = f'rsem-calculate-expression --paired-end {paired[0]} {paired[1]} --bowtie2 -p {threads} --no-bam-output {output}{prefix}/Trinity {output}{prefix}/{prefix}'
        #     # print(calculate_cmd)
        #     os.system(calculate_cmd)
        for paired in paireds:
            # print(paired)
            if not is_pipeline:
                templist = list(paired)
                for i in range(len(templist)):
                    filename = os.path.basename(templist[i])  # file.txt
                    templist[i] = filename
                prefix = longest_common_prefix(templist).strip("_|.")
            else:
                inverse_dict = {tuple(v): k for k, v in trim_files.items()}
                key = paired
                if key in inverse_dict:
                    prefix = inverse_dict[key]
                else:
                    key = (key[1],key[0])
                    prefix = inverse_dict[key]
            # os.makedirs(f"{output}{prefix}/{prefix}")
            # print(prefix)
            calculate_cmd = f'rsem-calculate-expression --paired-end {paired[0]} {paired[1]} --bowtie2 -p {threads} --no-bam-output {output}Trinity {output}{prefix}'
            # print(calculate_cmd)
            reslist.append(f"{output}{prefix}.genes.results")
            rc = os.system(calculate_cmd)
            # 判断是否执行成功
            if rc == 0:
                logger.info(f"rsem-calculate-expression命令执行成功！\n")
            else:
                logger.error(f"rsem-calculate-expression命令执行失败！\n")
                raise Exception(f"rsem-calculate-expression命令执行失败！")
    # 调用系统命令执行rsem-calculate-expression
    # os.system(calculate_cmd)
    # print(calculate_cmd)
    # 打印结束信息
    logger.info(f"完成运行trinity rsem!\n")
    logger.info(f"开始整合表达数据...\n")
    merge_cmd = f'rsem-generate-data-matrix {" ".join(reslist)} > {output}genes.expected_count.matrix'
    rc = os.system(merge_cmd)
    if rc == 0:
        logger.info(f"rsem-generate-data-matrix命令执行成功！\n")
    else:
        logger.error(f"rsem-generate-data-matrixn命令执行失败！\n")
        raise Exception(f"rsem-generate-data-matrix命令执行失败！")

    # 读取原始 pkl 文件中的列表数据
    config_pkl = f'{config_prefix}.pkl'
    with open(config_pkl, "rb") as f:
        data = pickle.load(f)

    # 在原始列表中添加新的字符串元素
    data["expression_matrix"] = f"{output}genes.expected_count.matrix"

    # 将原始列表和新的字符串一起写入新的 pkl 文件中
    with open(config_pkl, "wb") as f:
        pickle.dump(data, f)

# 如果是主程序，则调用命令行函数
if __name__ == '__main__':
    rsem_quantify()