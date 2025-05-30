import sys

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
from .GC123 import correlationStat

def correlation(correlationstat,figname):
    # 生成示例数据（替换为实际数据读取）
    data = pd.read_csv(correlationstat, sep='\t')

    # 数据预处理
    for col in data.columns:
        data[col] = pd.to_numeric(data[col], errors='coerce')

    # 计算相关系数和P值矩阵
    corr_matrix = data.corr(method='pearson').round(2)
    p_matrix = np.zeros_like(corr_matrix, dtype=float)

    for i in range(data.shape[1]):
        for j in range(data.shape[1]):
            if i != j:
                valid_idx = data.iloc[:, [i,j]].dropna().index
                x = data.iloc[valid_idx, i]
                y = data.iloc[valid_idx, j]
                _, p = stats.pearsonr(x, y)
                p_matrix[i, j] = p

    # 创建注释矩阵
    annot_matrix = np.empty_like(corr_matrix, dtype=object)
    for i in range(corr_matrix.shape[0]):
        for j in range(corr_matrix.shape[1]):
            if i == j:
                annot_matrix[i, j] = "1.00"  # 强制显示对角线为1.00
            elif i > j:
                annot_matrix[i, j] = f"{corr_matrix.iloc[i, j]:.2f}"  # 下半部分显示数值
            else:
                # 上半部分显示显著性星号
                p = p_matrix[i, j]
                if p < 0.001:
                    annot_matrix[i, j] = "***"
                elif p < 0.01:
                    annot_matrix[i, j] = "**"
                elif p < 0.05:
                    annot_matrix[i, j] = "*"
                else:
                    annot_matrix[i, j] = ""

    # 创建全矩阵mask（不隐藏任何单元格）
    mask = np.zeros_like(corr_matrix, dtype=bool)

    # 绘制热图
    plt.figure(figsize=(12, 10))
    sns.heatmap(corr_matrix,
                mask=mask,
                annot=annot_matrix,  # 使用自定义注释矩阵
                fmt='',              # 禁用自动格式化
                cmap='coolwarm',
                vmin=-1,
                vmax=1,
                center=0,
                square=True,
                linewidths=0.5,
                cbar_kws={"shrink": 0.8},
                annot_kws={"size": 12})

    # 调整坐标轴
    plt.xticks(rotation=45, ha='right', fontsize=12)
    plt.yticks(rotation=0, fontsize=12)
    plt.title("Enhanced Correlation Matrix with Full Color", pad=20, fontsize=14)

    # 保存PDF
    plt.savefig(figname,
                format='pdf',
                dpi=300,
                bbox_inches='tight')
def run(infasta,codonWout,figName,outputPath,figType = "pdf"):

    correlationstat = outputPath + "correlationstat.csv"
    correlationStat(infasta, correlationstat, codonWout)

    figname = figName + "." + figType
    correlation(correlationstat, figname)
if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python Correlation.py <inputfile> <outputfile> <figname>")
        sys.exit(1)
    correlationstat = "correlationstat.csv"
    infasta = sys.argv[1]
    codonWout = sys.argv[2]
    correlationStat(infasta, correlationstat, codonWout)
    figname = sys.argv[3]
    correlation(correlationstat,figname)