import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import Seq


def calculate_gc(seq):
    """手动计算GC含量百分比"""
    gc = seq.count('G') + seq.count('C')
    total = len(seq)
    return round((gc / total) * 100, 2) if total > 0 else 0.0



def calculate_GC123(seq):
    GC1 = []
    GC2 = []
    GC3 = []
    for i in range(1, len(seq) + 1):
        if i % 3 == 1:
            GC1.append(str(seq[i - 1]))
        elif i % 3 == 2:
            GC2.append(str(seq[i - 1]))
        else:
            GC3.append(str(seq[i - 1]))
    GC1_seq = Seq(''.join(GC1))
    GC2_seq = Seq(''.join(GC2))
    GC3_seq = Seq(''.join(GC3))

    return [round(calculate_gc(GC1_seq), 2), round(calculate_gc(GC2_seq), 2), round(calculate_gc(GC3_seq), 2)]

def correlationStat(in_fasta,out_results,codonWout):
    # 准备列名
    columns = ["title", "GC1", "GC2", "GC3", "GC"]

    # 创建空 DataFrame
    data = []

    # 遍历 FASTA 文件并收集数据
    for rec in SeqIO.parse(in_fasta, 'fasta'):
        # 计算 GC_total（保留两位小数）
        GC_total = round(calculate_gc(rec.seq), 2)

        # 计算 GC1, GC2, GC3
        GC1, GC2, GC3 = calculate_GC123(rec.seq)

        # 将数据添加到列表中（注意顺序要与 columns 对应）
        data.append([
            rec.id,  # gene_name
            str(GC1),  # GC1
            str(GC2),  # GC2
            str(GC3),  # GC3
            str(GC_total)  # GC_total
        ])

    # 转换为 DataFrame
    df1 = pd.DataFrame(data, columns=columns)
    df2 = pd.read_csv(codonWout, sep='\t', usecols=['title','CAI','CBI','Fop','Nc','Gravy','Aromo'])
    merged_df = pd.merge(df1, df2, on='title', how='inner')
    del merged_df['title']

    merged_df.to_csv(out_results, sep='\t', index=False)
