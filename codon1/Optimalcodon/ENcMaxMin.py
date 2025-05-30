import sys

import pandas as pd
import math
from Bio import SeqIO
from .RSCUofCDS import calculate_overall_rscu
def maxMinENc(inputfile,rat=0.1):
    '''
    inputfile = 'Psal.out'
    min_10_percent = 'Psal.min_10_percent.txt'
    max_10_percent = 'Psal.max_10_percent.txt'
    '''

    # 读取文件（假设为制表符分隔）
    df = pd.read_csv(inputfile, sep='\t')

    # 处理Nc列：转换为数值，非数值设为NaN
    df['Nc'] = pd.to_numeric(df['Nc'], errors='coerce')

    # 删除包含非数值的行
    df_clean = df.dropna(subset=['Nc']).copy()

    # 按Nc值排序
    sorted_df = df_clean.sort_values(by='Nc')

    # 计算10%的数量（向上取整）
    total = len(sorted_df)
    n = math.ceil(total * rat)

    # 提取最小和最大的10%
    min_10 = list(sorted_df.head(n)['title'])
    max_10 = list(sorted_df.tail(n)['title'])


    # 保存结果（不保留标题行）
    # min_10.to_csv(min_10_percent, index=False, header=False)
    # max_10.to_csv(max_10_percent, index=False, header=False)
    #
    print(f"处理完成：共{total}条有效数据，各取{n}条")
    return min_10, max_10

def extract_sequences(fasta_file, names_list, output_file):
    # # 读取序列名称列表
    # with open(names_file, 'r') as names_fh:
    #     names = [line.strip() for line in names_fh]
    names = names_list
    # 提取匹配的序列
    extracted_sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id in names:
            extracted_sequences.append(record)

    # 保存提取的序列到新的FASTA文件
    with open(output_file, 'w') as output_fh:
        SeqIO.write(extracted_sequences, output_fh, "fasta")

def optimalstat(cds,codonWout,outstat):
    min_10, max_10 = maxMinENc(codonWout)
    fasta_file = cds
    # 提取序列
    extract_sequences(fasta_file, min_10, 'min.fasta')
    extract_sequences(fasta_file, max_10, 'max.fasta')
    # 计算RSCU
    min_rscu_values = calculate_overall_rscu('min.fasta')
    max_rscu_values = calculate_overall_rscu('max.fasta')
    all_rscu_values = calculate_overall_rscu(fasta_file)
    # 转换为DataFrame输出
    dfmin = pd.DataFrame.from_dict(min_rscu_values, orient='index', columns=['min10RSCU'])
    dfmax = pd.DataFrame.from_dict(max_rscu_values, orient='index', columns=['max10RSCU'])

    dfall = pd.DataFrame.from_dict(all_rscu_values, orient='index', columns=['RSCU'])
    dfmin.sort_index(inplace=True)
    dfmax.sort_index(inplace=True)
    dfall.sort_index(inplace=True)
    dfmax.index.name = 'Codon'
    dfmin.index.name = 'Codon'
    dfall.index.name = 'Codon'
    merged_df = dfmin.join(dfmax).join(dfall)
    merged_df['D-value'] = merged_df['max10RSCU'] - merged_df['min10RSCU']
    merged_df.to_csv(outstat, sep='\t')
    return merged_df

if __name__ == '__main__':
    cds = "Psal.cds.fasta"
    codonWout = "Psal.out"
    outstat = "optimalCodon.csv"
    optimalstat(cds,codonWout,outstat)
    # min_10, max_10 = maxMinENc('Psal.out')
    # fasta_file = "Psal.cds.fasta"
    # # 提取序列
    # extract_sequences(fasta_file, min_10, 'min.fasta')
    # extract_sequences(fasta_file, max_10, 'max.fasta')
    # # 计算RSCU
    # min_rscu_values = calculate_overall_rscu('min.fasta')
    # max_rscu_values = calculate_overall_rscu('max.fasta')
    # all_rscu_values = calculate_overall_rscu(fasta_file)
    # # 转换为DataFrame输出
    # dfmin = pd.DataFrame.from_dict(min_rscu_values, orient='index', columns=['min10RSCU'])
    # dfmax = pd.DataFrame.from_dict(max_rscu_values, orient='index', columns=['max10RSCU'])
    #
    # dfall = pd.DataFrame.from_dict(all_rscu_values, orient='index', columns=['RSCU'])
    # dfmin.sort_index(inplace=True)
    # dfmax.sort_index(inplace=True)
    # dfall.sort_index(inplace=True)
    # dfmax.index.name = 'Codon'
    # dfmin.index.name = 'Codon'
    # dfall.index.name = 'Codon'
    # merged_df = dfmin.join(dfmax).join(dfall)
    # merged_df['D-value'] = merged_df['max10RSCU'] - merged_df['min10RSCU']
    # merged_df.to_csv("optimalCodon.csv")
