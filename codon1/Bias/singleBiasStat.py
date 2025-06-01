from Bio import SeqIO
from collections import defaultdict
import csv

def singleBiasStat(input_cds,outstat):

    # 目标氨基酸密码子字典（需全大写）
    AMINO_ACID_CODONS = {
        'Alanine': ['GCT', 'GCC', 'GCA', 'GCG'],
        'Arginine': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
        'Glycine': ['GGT', 'GGC', 'GGA', 'GGG'],
        'Leucine': ['CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'],
        'Proline': ['CCT', 'CCC', 'CCA', 'CCG'],
        'Serine': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
        'Threonine': ['ACT', 'ACC', 'ACA', 'ACG'],
        'Valine': ['GTT', 'GTC', 'GTA', 'GTG']
    }

    # 创建目标密码子集合
    target_codons = set()
    for codons in AMINO_ACID_CODONS.values():
        target_codons.update([c.upper() for c in codons])

    # 准备结果容器
    results = []

    # 遍历每条基因序列
    for record in SeqIO.parse(input_cds, "fasta"):
        gene_id = record.id
        seq = str(record.seq).upper()

        # 初始化计数器
        counter = defaultdict(int)
        total = 0

        # 检查序列长度有效性
        if len(seq) % 3 != 0:
            print(f"警告：基因 {gene_id} 长度不是3的倍数，已跳过")
            continue

        # 遍历所有密码子
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]

            # 检查密码子有效性
            if len(codon) != 3 or any(n not in "ATGC" for n in codon):
                continue

            # 判断是否为目标密码子
            if codon in target_codons:
                third_base = codon[2]
                counter[third_base] += 1
                total += 1

        # 计算百分比
        percentages = {}
        for base in ['A', 'T', 'G', 'C']:
            count = counter.get(base, 0)
            percentages[f"{base}%"] = round((count / total * 100), 2) if total > 0 else 0.0

        # 记录结果
        results.append({
            "GeneID": gene_id,
            "Total_Codons": total,
            **percentages
        })

    # 写入CSV文件
    with open(outstat, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["GeneID", "Total_Codons", "A%", "T%", "G%", "C%"])
        writer.writeheader()
        writer.writerows(results)
