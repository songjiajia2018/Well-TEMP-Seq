#!/usr/bin/env python
#coding=utf-8

import re
import codecs
import argparse

parser = argparse.ArgumentParser(description='count id_duplication in a fastq_file and remove them from th file: \
                                              python fastq_remove_id_duplication.py \
                                              -F1 fastq_file1 \
                                              -F2 fastq_file2')

parser.add_argument('-F1', '--FASTQ1', type=str, required=True, help='Input fastq file for single end data, or first read in paired end')
parser.add_argument('-F2', '--FASTQ2', type=str, default='', help='Input fastq file for the second read of paired end data. Default value: null.')

args = parser.parse_args()


def rm_dup_id(file_url):
    file=open(file_url,'r')

    ##extract id
    list=[]
    pattern=r'@[A-Z0-9:-]*?\s'
    for i in file:
        if (re.search(pattern,i)):
            m=re.search(pattern,i)
            id=m.group(0)
            list.append(id)
    
    ##id_duplication.txt
    txt=codecs.open('id_duplication.txt','w')
    count=0
    list2=[]
    dict={}
    for i in list:
        if i in dict:
            dict[i]=dict[i]+1
        else:
            dict[i]=1
    for k, v in dict.items():
        if (v > 1):
            count=count+1
            c=str(count)
            n=str(v)
            txt.write(c+'\t'+k+'\t'+n+'\n')
            list2.append(k)
    txt.close()

    file.close()

    ##remove id_duplication from fastq_file
    file=open(file_url,'r')
    cont=file.read()
    file_name=file_url.split('/')[-1].split('.')[0]
    file_format=file_url.split('/')[-1].split('.')[-1]
    fastq_rm=codecs.open(file_name+'_remove_id_duplication.'+file_format,'w')
    for i in list2:
        print(i)
        pattern=i+'.*?\n'+'.*?\n'+'.*?\n'+'.*?\n'
        cont=re.sub(pattern,'',cont)
    fastq_rm.write(cont)
    fastq_rm.close()
    file.close()

rm_dup_id(args.FASTQ1)
if (args.FASTQ2):
    rm_dup_id(args.FASTQ2)
