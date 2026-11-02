#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@File    :   Uniprot_query.py
@Contact :   2454888366@qq.com

@Modify Time      @Author    @Version    @Desciption
------------      -------    --------    -----------
2023/3/23 上午11:17   skychou      1.0         None
'''
import os
import pickle
import re
import time
import urllib.parse
import logging as lg
import click
# 导入所需的模块
import requests
###分页操作
from requests.adapters import HTTPAdapter, Retry

# https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=fasta&query=%28organism_id%3A2697049%29%20AND%20%28reviewed%3Atrue%29
# https://rest.uniprot.org/uniprotkb/search?query=reviewed:true+AND+organism_id:{organism_id}&format=tsv

re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def get_batch(batch_url):
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        total = response.headers["x-total-results"]
        yield response, total
        batch_url = get_next_link(response.headers)


# 定义一个函数，用于发送请求并获取响应
def get_uniprot_data(query, format, fields: list, file):
    global data_success
    # try:
    # https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=fasta&query=%28organism_id%3A2697049%29%20AND%20%28reviewed%3Atrue%29
    UniprotKB_API = "https://rest.uniprot.org/uniprotkb/search"
    # 发送请求并获取响应
    params = {
        "query": query,
        "format": format,
        "fields": ",".join(fields),
        "size": 500
    }
    encoded_params = urllib.parse.urlencode(params)
    url = UniprotKB_API + "?" + encoded_params
    # 发送请求并获取响应
    progress = 0
    print("正在保存summary文件......")
    with open(file, 'w') as f:
        for batch, total in get_batch(url):
            lines = batch.text.splitlines()
            if not progress:
                print(lines[0], file=f)
            for line in lines[1:]:
                print(line, file=f)
            progress += len(lines[1:])
            print(f'{progress} / {total}')
    # print(response.url)
    # 'query': '(organism_id:2697049) AND (reviewed:true)'
    # 检查响应状态码是否正常
    data_success = True
    # except:
    #     print("出错啦!")
    #     data_success = False


def get_uniprot_fasta(query, file):
    global fasta_success
    try:
        # https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=fasta&query=%28organism_id%3A2697049%29%20AND%20%28reviewed%3Atrue%29
        UniprotKB_API = "https://rest.uniprot.org/uniprotkb/stream"
        params = {
            "compressed": "true",
            "format": "fasta",
            "query": query,
        }
        encoded_params = urllib.parse.urlencode(params)
        url = UniprotKB_API + "?" + encoded_params
        print("正在保存fasta文件......")
        with requests.get(url, stream=True) as request:
            request.raise_for_status()
            with open(file, 'wb') as f:
                for chunk in request.iter_content(chunk_size=2 ** 20):
                    f.write(chunk)
        fasta_success = True
    except:
        fasta_success = False


# 定义一个函数，用于记录下载信息到本地文件
def record_download_info(file, query, file_prefix, addition_fields):
    # 打开本地文件
    with open(file, "a") as f:
        # 写入检索关键词和时间
        f.write(f">")
        f.write(f"\nQuery:{query}")
        f.write(f"\nFasta_file:{file_prefix}.fasta")
        f.write(f"\nSummary_file:{file_prefix}.summary")
        f.write(f"\nAddition_fields:{addition_fields}")
        f.write(f"\nTime:{time.ctime()}\n")


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


# 定义一个命令行接口，使用click模块
@click.command()
@click.option("--query", "-q",
              help='检索关键词，例如"(taxonomy_id:4751) AND (reviewed:true)"\n详细请查看手册：https://www.uniprot.org/help/query-fields',
              default="")
# @click.option("--size", "-s", default="100", help="检索该物种信息的数目(行)，默认为100条")
@click.option("--addition_fields", "-af", default="",
              help='检索该物种信息的额外条目(列)，用逗号分割，默认为:"accession", "xref_geneid","gene_names", "organism_id", "organism_name", "cc_pathway", "go_id", "go"')
@click.option("--format", "-f", default="tsv", help="指定下载总结文件的格式")
@click.option("--output_dir", "-o", default="./", help="结果输出目录，默认为当前目录")
@click.option("--fasta_output", "-fo", default="", help="fasta文件名称，默认为{query}.fasta")
@click.option("--summary_output", "-so", default="", help="summary表格文件名称，默认为{query}.tsv")
@click.option("--uniprotKB_records", "-u", default="UniprotKB_records.txt",
              help="下载记录文件，默认为UniprotKB_records.txt")
@click.option("--config", "-c", type=click.Path(), default="",
              help="配置文件目录")
def main(query, addition_fields, format, output_dir,
         fasta_output, summary_output, uniprotkb_records, config):
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

    # 检查输入参数是否有效
    if not query:
        print("Please provide query.")
        return
    if format not in ["tsv"]:
        print("Please provide a valid format.")
        return
    file_prefix = format_filename(query)
    if fasta_output == "":
        fasta_output = file_prefix + ".fasta.gz"
    if summary_output == "":
        summary_output = file_prefix + ".summary"

    # 获取data的准备工作
    default_field = ["accession", "xref_geneid", "gene_names", "lineage", "lineage_ids", "organism_name", "organism_id",
                     "go_id", "go", "go_p", "go_c", "go_f", "cc_pathway"]

    if addition_fields != "":
        try:
            addition_fields = addition_fields.split(",")
            default_field.extend(addition_fields)
        except Exception as e:
            raise Exception(f"请用,作为分割符！{e}")
    fields = default_field

    # 检查下载文件是否存在，判断是否进行后续下载步骤
    data_flag = True
    fasta_flag = True
    uniprotkb_records = os.path.join(output_dir, uniprotkb_records)
    if os.path.exists(uniprotkb_records):
        with open(uniprotkb_records, 'r') as f:
            text = f.read()
        records = text.split(">\n")
        check_addition_fields = f"Addition_fields:{addition_fields}"
        check_Query = f"Query:{query}"
        for i in records:
            if i != "" and check_Query in i:
                fasta_flag = False
                break
        for i in records:
            if not fasta_flag and check_addition_fields in i:
                data_flag = False
                break
        if os.path.exists(os.path.join(output_dir, summary_output)) and not data_flag:
            data_flag = False
        if os.path.exists(os.path.join(output_dir, fasta_output)) and not fasta_flag:
            fasta_flag = False

    if data_flag:
        # 调用函数保存summary文件
        try:
            get_uniprot_data(query, format, fields, os.path.join(output_dir, summary_output))
            print(f"summary结果保存为{os.path.join(output_dir, summary_output)}")
        except Exception as e:
            print(e)
    else:
        print(f"summary文件已经存在！为{os.path.join(output_dir, summary_output)}")

    if fasta_flag:
        # 调用函数保存fasta文件
        try:
            get_uniprot_fasta(query, os.path.join(output_dir, fasta_output))
            print(f"fasta结果保存为{os.path.join(output_dir, fasta_output)}")
        except Exception as e:
            print(e)
    else:
        print(f"fasta文件已经存在！为{os.path.join(output_dir, fasta_output)}")

    if data_flag and fasta_flag:
        flag = data_success and fasta_success
    elif data_flag:
        flag = data_success
    elif fasta_flag:
        flag = fasta_success
    else:
        flag = False
    if flag:
        record_download_info(uniprotkb_records, query, file_prefix, addition_fields)
        print("执行成功！")
    else:
        if not fasta_flag and not data_flag:
            print("执行结束，summary和fasta文件均存在！")
        else:
            raise Exception(f"执行失败！执行情况：是否获得fasta文件：{fasta_success}。是否获得summary文件：{data_success}。")

    if pipeline:
        # 读取原始 pkl 文件中的列表数据
        config_pkl = f'{config_prefix}.pkl'
        with open(config_pkl, "rb") as f:
            data = pickle.load(f)

        # 在原始列表中添加新的字符串元素
        tempdict = {"fasta_file": os.path.join(output_dir, fasta_output),
                    "summary_file": os.path.join(output_dir, summary_output)}
        data["uniport_query_files"] = tempdict

        # 将原始列表和新的字符串一起写入新的 pkl 文件中
        with open(config_pkl, "wb") as f:
            pickle.dump(data, f)

if __name__ == "__main__":
    main()
