import csv

def parse_fasta(file_path):
    """解析FASTA文件，返回基因ID与序列的字典"""
    genes = {}
    current_id = None
    current_seq = []

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id is not None:
                    genes[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]  # 取第一个字段作为ID
                current_seq = []
            else:
                current_seq.append(line.upper())

    if current_id is not None:
        genes[current_id] = ''.join(current_seq)

    return genes


def calculate_rscu(cds_sequence):
    # 定义标准遗传密码表
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

    # 构建氨基酸到密码子列表的映射
    amino_to_codons = {}
    for codon, amino in genetic_code.items():
        if amino == '*':
            continue  # 跳过终止密码子
        if amino not in amino_to_codons:
            amino_to_codons[amino] = []
        amino_to_codons[amino].append(codon)

    # 分割CDS序列为密码子并转换为大写
    codons = [cds_sequence[i:i + 3].upper() for i in range(0, len(cds_sequence), 3)]

    # 统计密码子出现次数（排除终止密码子和无效密码子）
    counts = {}
    stop_codons = {'TAA', 'TAG', 'TGA'}
    for codon in codons:
        if len(codon) != 3:
            continue  # 忽略长度不足的片段
        if codon in stop_codons:
            continue  # 过滤终止密码子
        amino = genetic_code.get(codon, None)
        if amino is None or amino == '*':
            continue  # 无效或终止密码子
        counts[codon] = counts.get(codon, 0) + 1

    # 计算每个密码子的RSCU值
    rscu = {}
    for codon in counts:
        amino = genetic_code[codon]
        family = amino_to_codons.get(amino, [])
        total = sum(counts.get(c, 0) for c in family)
        num_codons = len(family)
        if num_codons == 0 or total == 0:
            continue  # 避免除零错误（理论上不会发生）
        expected = total / num_codons
        rscu_val = counts[codon] / expected
        rscu[codon] = round(rscu_val, 4)  # 保留四位小数

    return rscu

def process_cds_file(input_file, output_file):
    """主处理函数"""
    # 1. 读取FASTA文件
    genes = parse_fasta(input_file)

    # 2. 收集所有可能的密码子作为CSV表头
    all_codons = set()
    results = {}
    for gene_id, seq in genes.items():
        rscu = calculate_rscu(seq)
        results[gene_id] = rscu
        all_codons.update(rscu.keys())

    # 3. 写入CSV文件
    header = ['GeneID'] + sorted(all_codons)
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(header)

        for gene_id, rscu in results.items():
            row = [gene_id]
            for codon in header[1:]:
                row.append(f"{rscu.get(codon, 0):.4f}")  # 无数据则填0
            writer.writerow(row)
