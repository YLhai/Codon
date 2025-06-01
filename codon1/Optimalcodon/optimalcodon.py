import pandas as pd
import matplotlib.pyplot as plt
from .ENcMaxMin import optimalstat



def optimalCodonPlot(optimalCodonStat,figname):

    # --- 数据加载 ---
    df = pd.read_csv(optimalCodonStat,sep='\t')

    # --- 筛选关键密码子 ---
    key_codons = df[(df['RSCU'] > 1) & (df['D-value'] > 0.08)]

    # --- 创建画布 ---
    plt.figure(figsize=(10, 8))
    ax = plt.gca()
    # --- 绘制全量数据点（灰色背景）---
    ax.scatter(
        x=df['D-value'],
        y=df['RSCU'],
        c='lightgrey',  # 非关键点用灰色
        alpha=0.5,
        s=40,
        label='Other Codons'
    )

    # --- 突出显示关键点 ---
    scatter = ax.scatter(
        x=key_codons['D-value'],
        y=key_codons['RSCU'],
        c='#FF6B6B',  # 关键点用红色
        edgecolors='w',
        linewidth=0.5,
        s=80,
        label='Optimal Candidates'
    )

    # --- 添加参考线 ---
    ax.axvline(0.08, color='dodgerblue', ls='--', lw=1.5, alpha=0.8)
    ax.axhline(1.0, color='mediumseagreen', ls='-.', lw=1.5, alpha=0.8)

    # --- 添加筛选标签 ---
    texts = []
    for idx, row in key_codons.iterrows():
        texts.append(ax.text(
            x=row['D-value'],
            y=row['RSCU'],
            s=row['Codon'],
            fontsize=10,
            weight='bold',  # 加粗标签
            color='#2D4263'
        ))

    # --- 格式美化 ---
    ax.set_xlabel('ΔRSCU (D-value)', fontsize=12)
    ax.set_ylabel('RSCU Value', fontsize=12)
    ax.set_title('Optimal Codon Candidates (RSCU>1 & ΔRSCU>0.08)', fontsize=14)
    ax.legend(loc='upper left', frameon=True)
    ax.grid(ls='--', alpha=0.3)

    # --- 保存输出 ---
    plt.tight_layout()
    plt.savefig(figname,
                format='pdf',
                dpi=300,
                bbox_inches='tight')
def run(cds,codonWout,figname,outpath,figType="pdf"):
    inputfile = outpath+"optimalCodon.csv"
    figname = figname + "." + figType
    optimalstat(cds,codonWout,inputfile)
    optimalCodonPlot(inputfile, figname)
