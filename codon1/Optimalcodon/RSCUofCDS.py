import pandas as pd
from collections import defaultdict


def calculate_overall_rscu(fasta_path):
    """计算包含终止密码子的整体RSCU值"""
    # 定义完整的遗传密码表（包含终止密码子）
    genetic_code = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'
    }

    # 构建氨基酸到密码子的映射（包含终止密码子）
    amino_to_codons = defaultdict(list)
    for codon, amino in genetic_code.items():
        amino_to_codons[amino].append(codon)

    # 统计所有密码子出现次数
    codon_counts = defaultdict(int)

    # 读取并处理FASTA文件
    with open(fasta_path) as f:
        current_seq = []
        for line in f:
            if line.startswith('>'):
                if current_seq:
                    process_sequence(''.join(current_seq), codon_counts, genetic_code)
                current_seq = []
            else:
                current_seq.append(line.strip().upper())
        if current_seq:  # 处理最后一个序列
            process_sequence(''.join(current_seq), codon_counts, genetic_code)

    # 计算RSCU值
    rscu = {}
    for codon, count in codon_counts.items():
        amino = genetic_code[codon]
        family = amino_to_codons[amino]
        total = sum(codon_counts[c] for c in family)
        if total == 0:
            continue
        expected = total / len(family)
        rscu[codon] = round(count / expected, 4) if expected != 0 else 0

    return rscu


def process_sequence(seq, counts, genetic_code):
    """处理单个序列的密码子计数"""
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        if len(codon) == 3 and codon in genetic_code:
            counts[codon] += 1

def run(fasta_file,RSCU_file):

    # 计算结果
    rscu_values = calculate_overall_rscu(fasta_file)
    # 转换为DataFrame输出
    df = pd.DataFrame.from_dict(rscu_values, orient='index', columns=['RSCU'])
    df.sort_index(inplace=True)
    df.index.name = 'Codon'
    df.to_csv(RSCU_file,sep='\t')
