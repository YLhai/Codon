import sys
from collections import defaultdict

# 定义目标氨基酸及其对应密码子
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

STOP_CODONS = {'TAA', 'TAG', 'TGA'}


def read_fasta(file_path):
    """读取FASTA文件生成器"""
    seq_id, seq = None, []
    with open(file_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq_id is not None:
                    yield (seq_id, ''.join(seq))
                seq_id = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.upper())
    if seq_id is not None:
        yield (seq_id, ''.join(seq))


def analyze_codons(input_file):
    """主分析函数"""
    # 初始化统计字典
    stats = defaultdict(
        lambda: {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'total': 0}
    )

    # 合并所有目标密码子
    target_codons = set()
    for codons in AMINO_ACID_CODONS.values():
        target_codons.update(codons)

    # 处理每条序列
    for seq_id, seq in read_fasta(input_file):
        # 分割密码子并过滤
        codons = [seq[i:i + 3] for i in range(0, len(seq), 3) if len(seq[i:i + 3]) == 3]

        for codon in codons:
            if codon in STOP_CODONS:
                continue

            if codon in target_codons:
                # 找到对应的氨基酸
                for aa, aa_codons in AMINO_ACID_CODONS.items():
                    if codon in aa_codons:
                        third_base = codon[2]
                        stats[aa][third_base] += 1
                        stats[aa]['total'] += 1
                        break

    return stats


def calculate_averages(stats):
    """计算八种氨基酸的平均值"""
    totals = {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'total': 0}
    for aa in AMINO_ACID_CODONS:
        if stats[aa]['total'] == 0:
            continue
        for base in ['A', 'T', 'G', 'C']:
            totals[base] += stats[aa][base]
        totals['total'] += stats[aa]['total']

    if totals['total'] == 0:
        return {base: 0.0 for base in ['A', 'T', 'G', 'C']}

    return {
        base: round((totals[base] / totals['total']) * 100, 2)
        for base in ['A', 'T', 'G', 'C']
    }


def write_results(stats, average, output_file):
    """写入结果文件"""
    with open(output_file, 'w') as f:
        # 写入表头
        f.write("AminoAcid\tA%\tT%\tG%\tC%\tTotalCodons\n")

        # 写入各氨基酸数据
        for aa in sorted(AMINO_ACID_CODONS.keys()):
            data = stats[aa]
            if data['total'] == 0:
                f.write(f"{aa}\t0.00\t0.00\t0.00\t0.00\t0\n")
                continue

            a_pct = round((data['A'] / data['total']) * 100, 2)
            t_pct = round((data['T'] / data['total']) * 100, 2)
            g_pct = round((data['G'] / data['total']) * 100, 2)
            c_pct = round((data['C'] / data['total']) * 100, 2)

            f.write(f"{aa}\t{a_pct}\t{t_pct}\t{g_pct}\t{c_pct}\t{data['total']}\n")

        # 写入平均值
        f.write("\nAverage\t")
        f.write(f"{average['A']}\t{average['T']}\t{average['G']}\t{average['C']}\t-\n")

def globBiasStat(input_cds, outStat):

    try:
        stats = analyze_codons(input_cds)
        average = calculate_averages(stats)
        write_results(stats, average, outStat)
        print(f"Analysis completed. Results saved to {outStat}")
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)
