#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@File    :   check.py    
@Contact :   2454888366@qq.com

@Modify Time      @Author    @Version    @Desciption
------------      -------    --------    -----------
10/9/2023 2:33 下午   skychou      1.0         None
'''
import sys
import subprocess
import os

if __name__ == '__main__':
    if sys.version_info < (3, 8):
        raise Exception("python({})版本过低,请更新至3.7以上版本！".format(sys.version))
    current_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(current_dir)
        # 尝试导入依赖
    try:
        import matplotlib
        matplotlib.use('Agg')
        import gseapy
        import numpy
        import glob
        import zipfile
        import shutil
        import logging as lg
        import os
        import shutil
        import pickle
        import click
        import rpy2
        import pandas
        import re
        import time
        import urllib
        import requests
    except ImportError as e:
        # 获取缺少的包名
        package = e.name
        # 调用pip命令来安装
        print(f"缺少python包:{package}，正在安装！如果安装失败请手动安装：pip install {package}")
        subprocess.run(["pip", "install", package])
        raise Exception(f"缺少python包:{package}，正在安装！如果安装失败请手动安装：pip install {package}")
    # 检查R以及R包
    try:
        from rpy2.robjects import r

        # 运行R的命令
        output = subprocess.check_output(['R', '--version'])
        version_string = output.decode('utf-8').split()[2]  # 从输出中获取版本字符串

        # 提取主版本号和次版本号
        major_version, minor_version = map(int, version_string.split('.')[:2])
    except Exception as e:
        raise Exception("rpy2导入失败，未能连接到R！")
    if major_version < 3 or (major_version == 3 and minor_version < 8):
        raise Exception("R({})版本过低,请更新至3.8以上版本！".format(major + "." + minor))

    from rpy2.robjects.packages import importr

    utils = importr('utils')

    def extract_versions(package_data):
        return dict(zip(
            package_data.rx(True, 'Package'),  # get Package column
            package_data.rx(True, 'Version')  # get Version column
        ))


    installed_packages = extract_versions(utils.installed_packages())
    r_packages = {'limma': '3.50.3', 'edgeR': '3.36.0',
                  'pheatmap': '1.0.12', 'data.table': '1.14.8'}  # 'DESeq2': '1.34.0',
    for packnames in r_packages.items():
        if not (packnames[0] in installed_packages):
            print(f"缺少R包:{packnames[0]}，正在安装！如果安装失败请手动安装：install.packages('{packnames[0]}')")
            utils.chooseCRANmirror(ind=1)  # select a CRAN mirror
            utils.install_packages(packnames[0])  # install the package
            raise Exception(
                f"缺少R包:{packnames[0]}，正在安装！如果安装失败请手动安装：install.packages('{packnames[0]}')")
        else:
            need_version = ".".join(packnames[1].split(".")[:-1])
            packages_version = ".".join(installed_packages[packnames[0]].split(".")[:-1])
            if float(packages_version) < float(need_version):
                print(
                    f"R包:{packnames[0]}版本过低，正在升级！如果升级失败请手动升级：library(remotes)install_version('{packnames[0]}', '{packnames[1]}')")
                utils.chooseCRANmirror(ind=1)  # select a CRAN mirror
                utils.install_packages(packnames[0])  # install the package
                raise Exception(
                    f"R包:{packnames[0]}版本过低，正在升级！如果升级失败请手动升级：library(remotes)install_version('{packnames[0]}', '{packnames[1]}')")

    # 检查和安装linux命令
    linux_packages = ["fastqc", "bowtie2", "pigz", "gzip"]
    for linux_package in linux_packages:
        code = subprocess.run(["which", linux_package]).returncode  # 检查fastqc命令是否存在
        if code != 0:
            print(f"缺少必备程序:{linux_package}，正在安装！如果安装失败请手动安装：apt-get install {linux_package}")
            raise Exception(f"缺少必备程序:{linux_package}!")
            # subprocess.run(["apt", "install", "fastqc"])  # 安装fastqc命令
            # subprocess.run(["apt", "update"])  # 更新软件源信息
            # subprocess.run(["apt", "upgrade"])  # 更新所有安装的ubuntu命令
    some_packages = ["TrimmomaticPE", "trimmomatic"]
    for linux_package in some_packages:
        code = subprocess.run(["which", linux_package]).returncode  # 检查fastqc命令是否存在
        if code == 0:
            break  #
    else:
        print(f"缺少必备程序:{linux_package}，正在安装！如果安装失败请手动安装：apt-get install {linux_package}")
        raise Exception(f"缺少必备程序:{linux_package}!")

    extend_software = ["Trinity", "RSEM"]
    for ext in extend_software:
        if ext == "Trinity":
            linux_package = "Trinity"
            code = subprocess.run(["which", linux_package]).returncode  # 检查fastqc命令是否存在
            if code != 0:
                if subprocess.run(["which", "conda"]).returncode:
                    subprocess.run(["conda", "install", "-c", "bioconda", "trinity"])
                    print(
                        f"缺少必备程序:{ext}，正在使用conda安装！如果安装失败请手动安装：conda install -c bioconda trinity")
                else:
                    print(f"缺少必备程序:{ext}，conda不存在！正在下载编译安装！如果安装失败请手动安装\n："
                          f"wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/Trinity-v2.15.1/trinityrnaseq-v2.15.1.FULL.tar.gz\n"
                          f"tar -xvzf trinityrnaseq-v2.15.1.FULL.tar.gz\n"
                          f"cd trinityrnaseq-v2.15.1\n"
                          f"make\n"
                          f"make install")
                    cmds = [
                        f"wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/Trinity-v2.15.1/trinityrnaseq-v2.15.1.FULL.tar.gz",
                        f"tar -xvzf trinityrnaseq-v2.15.1.FULL.tar.gz",
                        f"cd trinityrnaseq-v2.15.1",
                        f"make",
                        f"make install"
                    ]
                    for install_cmd in cmds:
                        try:
                            os.system(install_cmd)
                        except Exception as e:
                            raise Exception(f"安装{linux_package}失败！请手动安装！")
        if ext == "RSEM":
            linux_package = "rsem-prepare-reference"
            code = subprocess.run(["which", linux_package]).returncode  # 检查fastqc命令是否存在
            if code != 0:
                raise Exception(f"安装{ext}失败！请手动安装！")