#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@File    :   Quality_Control.py
@Contact :   2454888366@qq.com

@Modify Time      @Author    @Version    @Desciption
------------      -------    --------    -----------
3/7/2023 9:48 PM   skychou      1.0         None
'''

import glob
import logging as lg
import os
import pickle
import zipfile
import click

def run_fastqc(input_files, output_dir, format, threads):
    # 设置日志文件
    # 记录器
    global logger
    logger = lg.getLogger("mylog")
    logger.setLevel(lg.DEBUG)

    # 处理器
    consleHandler = lg.StreamHandler()
    consleHandler.setLevel(lg.DEBUG)

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

    # 没有给定handler日志级别，将用logger的级别
    global log_file
    log_file = f"{output_dir}quality_control.log"
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

    # 检查输入文件是否存在
    for file in input_files:
        if not os.path.exists(file):
            logger.error(f"输入文件{file}不存在，请检查路径")
            return

    # 检查输出目录是否存在，如果不存在则创建
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 转化为字符串
    input_files = " ".join(input_files)

    # 根据单端/双端参数设置fastqc命令
    fastqc_cmd = f"fastqc -o {output_dir} -f {format} --noextract --nogroup -t {threads} {input_files} >> {log_file}"
    flag = True  # 是否执行fastqc

    # 检测文件是否全部生成
    zip_files = glob.glob(f"{output_dir}*.zip")

    def get_pre1(filename):
        filename = os.path.basename(filename)  # file.txt
        prefix = os.path.splitext(filename)[0]  # file
        input_pre = prefix.strip("_fastqc")  # Con-1_1
        return input_pre

    def get_pre2(filename):
        filename = os.path.basename(filename)  # file.txt
        prefix = os.path.splitext(filename)[0]  # file
        suffix = os.path.splitext(filename)[1]
        if suffix == ".gz":
            prefix = os.path.splitext(prefix)[0]
        return prefix

    input_files_pre = [get_pre2(i) for i in input_files.split()]

    if zip_files != []:
        output_files_pre = [get_pre1(i) for i in zip_files]
        if set(input_files_pre) == set(output_files_pre):
            flag = False
    if flag:
        # 执行fastqc命令并打印结果
        logger.info(f"正在执行fastqc分析...\n")

        rc = os.system(fastqc_cmd)
        # 判断是否执行成功
        if rc == 0:
            logger.info(f"质量控制成功！\n")
        else:
            logger.error(f"质量控制失败！\n")
            raise Exception(f"质量控制失败！")
    else:
        logger.info(f"已经存在fastqc结果文件，将直接进行后续步骤！（如需要重新执行，请删除结果目录）！\n")


# 获得某些全局变量
def get_data(data_file):
    global bases #判断读长
    global gc_percent #GC百分比
    global all_base #所有碱基
    with zipfile.ZipFile(zip_file) as z:
        with z.open(data_file) as f:
            lines = f.read().decode("utf-8")
            lines = lines.split("\n")
            lines = [i for i in lines if i != ""]

    # 遍历文件中的每一行
    for line in lines:
        # 获得所有碱基数目
        if line.startswith("Total Sequences"):
            fields = line.split("\t")
            all_base = int(fields[1].strip())

        # 获得读长
        if line.startswith("Sequence length"):
            fields = line.split("\t")
            bases = int(fields[1].strip())

        # 获得GC含量
        if line.startswith("%GC"):
            fields = line.split("\t")
            gc_percent = int(fields[1].strip()) / 100



# 定义一个函数，用于计算Q20和Q30的百分比
def calculate_Q20_Q30(data_file):
    line_count = 0  # 记录行数
    per_sequence_quality = False  # 表示是否在每条序列质量值的部分
    total_q20 = 0  # 存储质量>=20的碱基数
    total_q30 = 0  # 存储质量>=30的碱基数
    with zipfile.ZipFile(zip_file) as z:
        with z.open(data_file) as f:
            lines = f.read().decode("utf-8")
            lines = lines.split("\n")
            lines = [i for i in lines if i != ""]

    # 遍历文件中的每一行
    for line in lines:
        # 增加行数计数器
        line_count += 1

        # 检查是否在每条序列质量值的部分
        if line.startswith(">>Per sequence quality scores"):
            per_sequence_quality = True

        elif line.startswith(">>END_MODULE"):
            per_sequence_quality = False

        elif per_sequence_quality:
            # 如果到了结尾
            if not per_sequence_quality:
                break
                # print((Q20_percent, Q30_percent))
                # 按制表符分割行，并获取质量值字符串（第二列）
            fields = line.split("\t")

            # 检查是否为表头
            if fields[0] == '#Quality':
                continue

            quality = float(fields[0].strip())
            count = float(fields[1].strip())

            # 计算total_q20和total_q30
            if quality >= 20:
                total_q20 += count
            if quality >= 30:
                total_q30 += count

    Q20_percent = round(total_q20 / all_base, 5)
    Q30_percent = round(total_q30 / all_base, 5)

    # 返回百分比的元组
    return (Q20_percent, Q30_percent)


# 定义一个函数，用于判断是否通过质量检测，规则是只有当三个指标全部为PASS或者WORN时通过检验。
def check_quality_status(summary_file):
    # 初始化一个空列表，用于存储三个指标的状态
    status_list = []

    # 打开summary文件，并按行读取内容
    with zipfile.ZipFile(zip_file) as z:
        with z.open(summary_file) as f:
            lines = f.read().decode("utf-8")
            lines = lines.split("\n")
            lines = [i for i in lines if i != ""]

    # 遍历每一行内容
    for line in lines:
        # 去掉换行符和空格，并按制表符分割字符串，得到一个列表
        line = line.strip().split('\t')

        # 判断是否是关于quality的指标，并将其状态添加到列表中（第一列为状态，第二列为指标）
        if line[1].startswith('Per base sequence quality') or line[1].startswith('Per tile sequence quality') or line[
            1].startswith('Per sequence quality scores'):
            status_list.append(line[0])

    # 判断列表中是否都是PASS或者WARN，如果是，则返回True，否则返回False。
    if all(status == 'PASS' or status == 'WARN' for status in status_list):
        return True

    else:
        return False


# 定义一个函数，用于判断当碱基ATCG含量趋于稳定时候的跳过的碱基数。
def find_stable_base_number(data_file, threshold):
    # 初始化一个空列表，用于存储每个位置的碱基含量
    base_content_list = []

    # 打开data文件，并按行读取内容
    with zipfile.ZipFile(zip_file) as z:
        with z.open(data_file) as f:
            lines = f.read().decode("utf-8")
            lines = lines.split("\n")
            lines = [i for i in lines if i != ""]

    # 遍历每一行内容
    for line in lines:
        # 去掉换行符和空格，并按制表符分割字符串，得到一个列表
        line = line.strip().split('\t')

        # 判断是否是关于Per base sequence content的数据，并将其添加到列表中（第一列为位置，后四列为A,T,C,G的百分比）
        if line[0].startswith('>>Per base sequence content') or (line[0].isdigit() and len(line) == 5):
            base_content_list.append(line)

    # 去掉第一行和最后一行，只保留数据部分
    base_content_list = base_content_list[1:-1]
    base_content_list_minus = base_content_list[::-1]

    def inter_fun(l):
        # 初始化一个变量，用于存储跳过的碱基数
        skip_number = 0

        # 遍历每个位置的数据
        for i in range(len(l)):
            # 将每个位置的A,T,C,G的百分比转换为浮点数，并计算其标准差（std）
            G = float(l[i][1])
            A = float(l[i][2])
            T = float(l[i][3])
            C = float(l[i][4])

            CG_std = (((C + G) / 100 - gc_percent) ** 2) ** 0.5

            # CG差(百分比)
            CG_m = ((C - G) ** 2) ** 0.5 / 100
            # AT差（百分比）
            AT_m = ((A - T) ** 2) ** 0.5 / 100

            # 判断三个标准差的值是否小于某个阈值（例如5），如果是，则认为该位置开始碱基含量趋于稳定，跳出循环；否则，累加跳过的碱基数。
            if CG_std < threshold and CG_m < threshold and AT_m < threshold:
                break

            else:
                skip_number += 1

        # 返回跳过的碱基数
        return skip_number

    skip_number = inter_fun(base_content_list)
    skip_number_minus = inter_fun(base_content_list_minus)

    return (skip_number, skip_number_minus)


# 定义一个总分析函数函数，用于调用上述三个函数，并打印结果。
def fastqc_parse(zip_file_name, out_dir, threshold):
    # 获得文件名
    filename = os.path.basename(zip_file_name)  # file.txt
    prefix = os.path.splitext(filename)[0]  # file

    res_lines = [f"file_name:{filename}\n"]
    # 全局变量，后面要用到
    global zip_file
    zip_file = zip_file_name
    # 使用zipfile模块获取其中的data和summary文件名
    with zipfile.ZipFile(zip_file_name) as z:
        data_file_name = z.namelist()[0] + 'fastqc_data.txt'
        summary_file_name = z.namelist()[0] + 'summary.txt'

    # 先获取所需全局变量
    get_data(data_file_name)

    # 调用calculate_Q20_Q30函数，传入data文件的名称，得到Q20和Q30的百分比
    Q20_percent, Q30_percent = calculate_Q20_Q30(data_file_name)

    # 打印Q20和Q30的百分比
    logger.info('Q20: {}'.format(Q20_percent))
    logger.info('Q30: {}'.format(Q30_percent))
    res_lines.append('Q20: {}\n'.format(Q20_percent))
    res_lines.append('Q30: {}\n'.format(Q30_percent))
    # 调用check_quality_status函数，传入summary文件的名称，得到是否通过质量检测的布尔值
    quality_status = check_quality_status(summary_file_name)

    # 打印是否通过质量检测的结果
    if quality_status:
        logger.info('通过质量检测！')
        res_lines.append('通过质量检测!True\n')

    else:
        logger.info('未通过质量检测！')
        res_lines.append('未通过质量检测!False\n')

    # 调用find_stable_base_number函数，传入data文件的名称，得到跳过的碱基数
    skip_numbers = find_stable_base_number(data_file_name, threshold)

    # 打印跳过的碱基数
    logger.info('开头减去的碱基数: {}'.format(skip_numbers[0]))
    res_lines.append('开头减去的碱基数: {}\n'.format(skip_numbers[0]))
    logger.info('结尾减去的碱基数: {}\n'.format(skip_numbers[1]))
    res_lines.append('结尾减去的碱基数: {}\n'.format(skip_numbers[1]))

    # 写到这里了，准备写保存文件和返回元组(quality_status,skip_number)!
    with open(f"{out_dir}{prefix}_summary.txt", "w+") as f:
        f.writelines(res_lines)

    return (quality_status, skip_numbers)


def trim(input_files: list, trimdict: dict, out_dir, software):
    trim_files = {}
    for file, trimtuple in trimdict.items():
        if not pairs:
            # 获得文件名
            filename = os.path.basename(file)  # file.txt
            prefix = os.path.splitext(filename)[0]  # file
            suffix = os.path.splitext(filename)[1]
            input_pre = prefix.strip("_fastqc")  # Con-1_1
            for i in input_files:
                if input_pre in i:
                    input_file = i
                    break


            output_file = os.path.join(out_dir, f"{input_pre}_trim.fq.gz")
        else:
            file1 = file[0]
            # 获得文件名
            filename = os.path.basename(file1)  # file.txt
            prefix = os.path.splitext(filename)[0]  # file
            suffix = os.path.splitext(filename)[1]
            input_pre = prefix.strip("_fastqc")  # Con-1_1
            for i in input_files:
                if input_pre in i:
                    input_file1 = i
                    break
            output_file1 = os.path.join(out_dir, f"{input_pre}_trim.fq.gz")
            other_output_file1 = os.path.join(out_dir, f"{input_pre}_trim_unpaired.fq.gz")
            input_pre1 = input_pre

            file2 = file[1]
            # 获得文件名
            filename = os.path.basename(file2)  # file.txt
            prefix = os.path.splitext(filename)[0]  # file
            suffix = os.path.splitext(filename)[1]
            input_pre = prefix.strip("_fastqc")  # Con-1_1
            for i in input_files:
                if input_pre in i:
                    input_file2 = i
                    break
            output_file2 = os.path.join(out_dir, f"{input_pre}_trim.fq.gz")
            other_output_file2 = os.path.join(out_dir, f"{input_pre}_trim_unpaired.fq.gz")
            input_pre2 = input_pre

        # 创建命令行
        # "cutadapt -j 0 -u 10 -u -10 -o output.fastq.gz input.fastq.gz"

        # 构造裁剪后的数据分组
        inverse_dict = {tuple(v): k for k, v in cfg.user_args["samples"].items()}
        if not pairs:
            key = (input_file,)
            trim_files[inverse_dict[key]] = [output_file]
        else:
            key = (input_file1, input_file2)
            trim_files[inverse_dict[key]] = [output_file1, output_file2]

        if software == "cutadapt":
            if not pairs:
                cmd = f"cutadapt -j 0 -u {trimtuple[1][0]} -u -{trimtuple[1][1]} -o {output_file} {input_file} "
            else:
                cmd = f"cutadapt -j 0 -u {trimtuple[0][1][0]} -u -{trimtuple[0][1][1]} -o {output_file1} {input_file1} && cutadapt -j 0 -u {trimtuple[0][1][0]} -u -{trimtuple[0][1][1]} -o {output_file2} {input_file2}"

        # 用到bases
        if software == "trimmomatic":

            if not pairs:
                if os.system("which trimmomatic")==0:
                    cmd = f"trimmomatic SE -threads {g_threads} {input_file} {output_file} HEADCROP:{trimtuple[1][0]} CROP:{int(bases) - int(trimtuple[1][0]) - int(trimtuple[1][1])}"
                else:
                    cmd = f"TrimmomaticSE -threads {g_threads} {input_file} {output_file} HEADCROP:{trimtuple[1][0]} CROP:{int(bases) - int(trimtuple[1][0]) - int(trimtuple[1][1])}"
            else:
                # "trimmomatic PE -threads 10 4m-1-LCM10095_combined_R1.fastq.gz 4m-1-LCM10095_combined_R2.fastq.gz trim_out/4m-1-LCM10095_combined_R1_trim.fastq.gz trim_out/4m-1-LCM10095_combined_R1_trim_unpaired.fastq.gz trim_out/4m-1-LCM10095_combined_R2_trim.fastq.gz trim_out/4m-1-LCM10095_combined_R2_trim_unpaired.fastq.gz HEADCROP:10 CROP:90"
                # cmd = f"trimmomatic PE -threads {g_threads} {input_file1} {input_file2} {output_file1} {output_file2} {other_output_file1} {other_output_file2} HEADCROP:{trimtuple[0][1][0]} CROP:{int(bases) - int(trimtuple[0][1][0]) - int(trimtuple[0][1][1])}"
                if os.system("which trimmomatic") == 0:
                    cmd = f"trimmomatic PE -threads {g_threads} {input_file1} {input_file2} {output_file1} {other_output_file1} {output_file2} {other_output_file2} HEADCROP:{trimtuple[0][1][0]} CROP:{int(bases) - int(trimtuple[0][1][0]) - int(trimtuple[0][1][1])}"
                else:
                    cmd = f"TrimmomaticPE -threads {g_threads} {input_file1} {input_file2} {output_file1} {other_output_file1} {output_file2} {other_output_file2} HEADCROP:{trimtuple[0][1][0]} CROP:{int(bases) - int(trimtuple[0][1][0]) - int(trimtuple[0][1][1])}"
        flag = True  # 是否执行剪切命令

        if not pairs:
            if os.path.exists(output_file):
                flag = False
            if flag:
                rc = os.system(cmd)
                # 判断是否执行成功
                if rc == 0:
                    logger.info(f"{input_pre}裁剪成功！裁剪长度为：{trimtuple[1][0]}，-{trimtuple[1][1]}\n")
                else:
                    logger.error(f"{input_pre}裁剪失败！\n")
                    raise Exception(f"{input_pre}裁剪失败！")
            else:
                logger.info(f"已经经过剪切！（如需要重新执行，请删除{output_file}）！\n")
        else:
            # print(cmd)
            if os.path.exists(output_file1) and os.path.exists(output_file2):
                flag = False
            if flag:
                rc = os.system(cmd)
                # 判断是否执行成功
                if rc == 0:
                    logger.info(
                        f"{input_pre1},{input_pre2}裁剪成功！裁剪长度为：{trimtuple[0][1][0]}，-{trimtuple[0][1][1]}\n")
                else:
                    logger.error(f"{input_pre1},{input_pre2}裁剪失败！\n")
                    raise Exception(f"{input_pre1},{input_pre2}裁剪失败！")
            else:
                logger.info(f"已经经过剪切！（如需要重新执行，请删除{output_file1,output_file2}）！\n")
    with open(f'{config_prefix}.pkl', 'wb') as f:
        data = {"trim_files":trim_files}
        pickle.dump(data, f)

@click.command()
@click.option("--input_files", "-i",default="",
              help="输入测序数据文件，多个文件用百分号(%)隔开。（双端测序文件左右端用逗号隔开）", type=click.STRING)
@click.option("--output_dir", "-o", default="./qc_fastqc_output", help="输出质控分析结果目录，默认为'fastqc_output'")
@click.option("--format", "-f", type=click.Choice(["fastq", "sam", "bam"]), default="fastq",
              help="测序数据文件格式,默认为'fastq'")
# @click.option("--paired", "-p", is_flag=True, help="是否为双端测序数据，默认为否")
@click.option("--not_cut", "-nc", is_flag=True, default=False, help="是否执行裁剪测序数据，默认为是")
@click.option("--threads", "-t", default="10", help="线程数，默认为10")
@click.option("--threshold", "-th", default=0.01, help="进行裁剪长度判断时的阈值，默认为0.01", type=float)
@click.option("--software", "-s", type=click.Choice(["trimmomatic", "cutadapt"]), default="trimmomatic",
              help="进行裁剪的程序，默认为trimmomatic")
@click.option("--config", "-c", type=click.Path(), default="./config.py",
              help="配置文件目录，默认为./config.py")
def main(input_files, output_dir, format, not_cut, threads, threshold, software, config):
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

    global g_threads
    g_threads = threads

    if not output_dir.endswith("/"):
        output_dir = output_dir + "/"
    if input_files=="":
        templist = []
        for files in cfg.user_args["samples"].values():
            templist.append(f"{files[0]},{files[1]}")
        input_files = "%".join(templist)

    # 将输入文件字符串转换为列表
    input_files = input_files.split("%")

    global pairs
    try:
        if "," in input_files[0]:
            pairs = True
        else:
            pairs = False
    except:
        logger.error(f"输入文件格式错误！{input_files}\n")
        raise Exception(f"输入文件格式错误！{input_files}")

    # 调用run_fastqc函数
    if not pairs:
        run_fastqc(input_files, output_dir, format, threads)
    else:
        realinput = []
        for i in input_files:
            realinput.extend(i.split(","))
        run_fastqc(realinput, output_dir, format, threads)

    # 调用fastqc_parse函数
    # 获得zip文件路径
    zip_files = glob.glob(f"{output_dir}*.zip")
    # zip_files: ['./qc_fastqc_output/Con-1_2_fastqc.zip', './qc_fastqc_output/Sam-2_1_fastqc.zip',...]
    trim_dict = {}
    for i in zip_files:
        filename = os.path.basename(i)  # file.txt
        prefix = os.path.splitext(filename)[0]  # file
        logger.info(f"{prefix} summary：")
        temp_res = fastqc_parse(i, output_dir, threshold)
        trim_dict[i] = temp_res  # './Con-1_2_fastqc.zip':(True,(10,10))
    if pairs:
        # 构建新的trim_dict，裁剪数值取最大
        new_trim_dict = {}
        for i in input_files:
            templist = i.split(",")  # ['./test/Con-1_1.fq.gz', './test/Con-1_2.fq.gz']
            for n in range(len(templist)):
                filename = templist[n]
                filename = os.path.basename(filename)  # file.txt
                prefix = os.path.splitext(filename)[0]  # file
                suffix = os.path.splitext(filename)[1]
                if suffix == ".gz":
                    prefix = os.path.splitext(prefix)[0]
                templist[n] = prefix
            keys = []
            for prefix in templist:
                for key in trim_dict.keys():
                    if prefix in key:
                        keys.append(key)
            max_plus = max(trim_dict[keys[0]][1][0], trim_dict[keys[1]][1][0])
            max_minus = max(trim_dict[keys[0]][1][1], trim_dict[keys[1]][1][1])
            new_trim_dict[(keys[0], keys[1])] = (trim_dict[keys[0]][0], (max_plus, max_minus)), (
            trim_dict[keys[1]][0], (max_plus, max_minus))

    if not not_cut:
        if not pairs:
            trim(input_files, trim_dict, output_dir, software)
        else:
            trim(realinput, new_trim_dict, output_dir, software)
    logger.info(f"Quality_Control完成！\n")


if __name__ == "__main__":
    main()
    # "python fastq_qc.py -i sample1.fastq,sample2.fastq -o qc_result"
    # 这个命令会对sample1.fastq和sample2.fastq两个双端测序数据进行质控分析，并且将结果保存在qc_result目录中，并且无视质控分析结果，执行后续步骤。
