import os
import pandas as pd
import matplotlib.pyplot as plt
from ..Optimalcodon import RSCUofCDS


def run(cds,figname,outputdir,figType = "pdf"):
    rscufile = os.path.join(outputdir,'RSCU.csv')
    figname = figname+"."+figType
    RSCUofCDS.run(cds,rscufile)
    rscu(rscufile,figname)

def rscu(rscu,figname):

    # 读取数据
    rscu_data = pd.read_csv(rscu,sep='\t')
    # codon_data = pd.read_csv(colorSet)
    # 硬编码颜色数据
    codon_data = pd.DataFrame([
    # 完整颜色编码数据
    {'Codon': 'TTT', 'AA': 'Phe', 'color': 6},
    {'Codon': 'TTC', 'AA': 'Phe', 'color': 5.5},
    {'Codon': 'TTA', 'AA': 'Leu', 'color': 6},
    {'Codon': 'TTG', 'AA': 'Leu', 'color': 5.5},
    {'Codon': 'TCT', 'AA': 'Ser', 'color': 6},
    {'Codon': 'TCC', 'AA': 'Ser', 'color': 5.5},
    {'Codon': 'TCA', 'AA': 'Ser', 'color': 5},
    {'Codon': 'TCG', 'AA': 'Ser', 'color': 4.5},
    {'Codon': 'TAT', 'AA': 'Tyr', 'color': 6},
    {'Codon': 'TAC', 'AA': 'Tyr', 'color': 5.5},
    {'Codon': 'TGT', 'AA': 'Cys', 'color': 6},
    {'Codon': 'TGC', 'AA': 'Cys', 'color': 5.5},
    {'Codon': 'TGG', 'AA': 'Trp', 'color': 6},
    {'Codon': 'CTT', 'AA': 'Leu', 'color': 5},
    {'Codon': 'CTC', 'AA': 'Leu', 'color': 4.5},
    {'Codon': 'CTA', 'AA': 'Leu', 'color': 4},
    {'Codon': 'CTG', 'AA': 'Leu', 'color': 3.5},
    {'Codon': 'CCT', 'AA': 'Pro', 'color': 6},
    {'Codon': 'CCC', 'AA': 'Pro', 'color': 5.5},
    {'Codon': 'CCA', 'AA': 'Pro', 'color': 5},
    {'Codon': 'CCG', 'AA': 'Pro', 'color': 4.5},
    {'Codon': 'CAT', 'AA': 'His', 'color': 6},  # 修正HIS为His
    {'Codon': 'CAC', 'AA': 'His', 'color': 5.5},
    {'Codon': 'CAA', 'AA': 'Gln', 'color': 6},
    {'Codon': 'CAG', 'AA': 'Gln', 'color': 5.5},
    {'Codon': 'CGT', 'AA': 'Arg', 'color': 6},
    {'Codon': 'CGC', 'AA': 'Arg', 'color': 5.5},
    {'Codon': 'CGA', 'AA': 'Arg', 'color': 5},
    {'Codon': 'CGG', 'AA': 'Arg', 'color': 4.5},
    {'Codon': 'ATT', 'AA': 'Ile', 'color': 6},
    {'Codon': 'ATC', 'AA': 'Ile', 'color': 5.5},
    {'Codon': 'ATA', 'AA': 'Ile', 'color': 5},
    {'Codon': 'ATG', 'AA': 'Met', 'color': 6},
    {'Codon': 'ACT', 'AA': 'Thr', 'color': 6},
    {'Codon': 'ACC', 'AA': 'Thr', 'color': 5.5},
    {'Codon': 'ACA', 'AA': 'Thr', 'color': 5},
    {'Codon': 'ACG', 'AA': 'Thr', 'color': 4.5},
    {'Codon': 'AAT', 'AA': 'Asn', 'color': 6},
    {'Codon': 'AAC', 'AA': 'Asn', 'color': 5.5},
    {'Codon': 'AAA', 'AA': 'Lys', 'color': 6},
    {'Codon': 'AAG', 'AA': 'Lys', 'color': 5.5},
    {'Codon': 'AGT', 'AA': 'Ser', 'color': 4},
    {'Codon': 'AGC', 'AA': 'Ser', 'color': 3.5},
    {'Codon': 'AGA', 'AA': 'Arg', 'color': 4},
    {'Codon': 'AGG', 'AA': 'Arg', 'color': 3.5},
    {'Codon': 'GTT', 'AA': 'Val', 'color': 6},
    {'Codon': 'GTC', 'AA': 'Val', 'color': 5.5},
    {'Codon': 'GTA', 'AA': 'Val', 'color': 5},
    {'Codon': 'GTG', 'AA': 'Val', 'color': 4.5},
    {'Codon': 'GCT', 'AA': 'Ala', 'color': 6},
    {'Codon': 'GCC', 'AA': 'Ala', 'color': 5.5},
    {'Codon': 'GCA', 'AA': 'Ala', 'color': 5},
    {'Codon': 'GCG', 'AA': 'Ala', 'color': 4.5},
    {'Codon': 'GAT', 'AA': 'Asp', 'color': 6},
    {'Codon': 'GAC', 'AA': 'Asp', 'color': 5.5},
    {'Codon': 'GAA', 'AA': 'Glu', 'color': 6},
    {'Codon': 'GAG', 'AA': 'Glu', 'color': 5.5},
    {'Codon': 'GGT', 'AA': 'Gly', 'color': 6},
    {'Codon': 'GGC', 'AA': 'Gly', 'color': 5.5},
    {'Codon': 'GGA', 'AA': 'Gly', 'color': 5},
    {'Codon': 'GGG', 'AA': 'Gly', 'color': 4.5},
    {'Codon': 'TAA', 'AA': 'Ter', 'color': 6},
    {'Codon': 'TAG', 'AA': 'Ter', 'color': 5.5},
    {'Codon': 'TGA', 'AA': 'Ter', 'color': 5}
])

    # 合并数据
    combined = pd.merge(rscu_data, codon_data, on="Codon")

    # 数据预处理
    combined['AA'] = pd.Categorical(combined['AA'])
    combined = combined.sort_values(['AA', 'color'])

    # 创建图形
    fig = plt.figure(figsize=(12, 8))
    gs = fig.add_gridspec(2, 1, height_ratios=[4, 1], hspace=0.05)

    # --------------------------
    # 主图 (堆叠柱状图)
    # --------------------------
    ax1 = fig.add_subplot(gs[0])

    # 按氨基酸分组（处理分类数据）
    groups = combined.groupby('AA', observed=True)
    aa_categories = combined['AA'].cat.categories  # 获取分类顺序

    # 初始化底部位置（使用浮点类型）
    bottom = pd.Series(0.0, index=aa_categories)

    # 颜色映射
    color_palette = {
        '6.0': '#8ea3c2',
        '5.5': '#a3b38c',
        '5.0': '#edb17f',
        '4.5':'#c9c9c9',
        '4.0': '#f1df82',
        '3.5': '#efafb4'
    }
    colors = {str(c): color_palette.get(str(c), '#999999') for c in combined['color'].unique()}

    # # 绘制堆叠条形（使用数值索引）
    # for idx, aa in enumerate(aa_categories):
    #     group = groups.get_group(aa)
    #     color_groups = group.groupby('color', observed=True)
    #     for color, sub_group in color_groups:
    #         ax1.bar(idx, sub_group['RSCU'].sum(),
    #                 bottom=bottom[aa],
    #                 color=colors[str(color)],
    #                 width=0.8)
    #         bottom[aa] += sub_group['RSCU'].sum()
    # 绘制堆叠条形（使用数值索引）
    for idx, aa in enumerate(aa_categories):
        group = groups.get_group(aa)

        # 按color数值降序排列（新添加的排序逻辑）
        sorted_colors = sorted(
            group['color'].unique(),
            key=lambda x: float(x),
            reverse=True
        )

        # 按排序后的颜色顺序绘制
        for color in sorted_colors:
            sub_group = group[group['color'] == color]
            ax1.bar(idx, sub_group['RSCU'].sum(),
                    bottom=bottom[aa],
                    color=colors[str(color)],
                    width=0.8)
            bottom[aa] += sub_group['RSCU'].sum()

    # 设置分类轴标签
    ax1.set_xticks(range(len(aa_categories)))
    ax1.set_xticklabels(aa_categories)
    ax1.tick_params(axis='x', which='both', length=0)
    ax1.set_ylabel('RSCU')
    ax1.set_ylim(0, 6.2)

    # --------------------------
    # 标签图
    # --------------------------
    # ax2 = fig.add_subplot(gs[1], sharex=ax1)
    #
    # # 创建标签位置（使用相同数值索引）
    # for idx, aa in enumerate(aa_categories):
    #     group = groups.get_group(aa)
    #     y_pos = 0.5
    #     for _, row in group.iterrows():
    #         ax2.text(idx, y_pos, row['Codon'],
    #                  ha='center', va='center',
    #                  fontsize=8,
    #                  bbox=dict(facecolor=colors[str(row['color'])],
    #                            edgecolor='none',
    #                            boxstyle='round,pad=0.2'))
    #         y_pos -= 0.2  # 调整标签垂直间距
    #
    # ax2.set_ylim(-1, 1)
    # ax2.axis('off')
    # 标签图
    ax2 = fig.add_subplot(gs[1], sharex=ax1)

    # 创建标签位置（与堆叠顺序一致）
    for idx, aa in enumerate(aa_categories):
        group = groups.get_group(aa)

        # 使用相同的颜色排序逻辑
        sorted_colors = sorted(
            group['color'].unique(),
            key=lambda x: float(x),
            reverse=True  # 保持与堆叠顺序一致
        )

        y_pos = 0.8  # 起始位置调高
        for color in sorted_colors:
            sub_group = group[group['color'] == color]
            for _, row in sub_group.iterrows():
                ax2.text(
                    idx, y_pos, row['Codon'],
                    ha='center', va='center',
                    fontsize=8,
                    bbox=dict(
                        facecolor=colors[str(color)],
                        edgecolor='none',
                        boxstyle='round,pad=0.2'
                    )
                )
                y_pos -= 0.25  # 向下排列

    ax2.set_ylim(-0.5, 1.2)
    ax2.axis('off')

    # 保存图形
    plt.tight_layout()
    plt.savefig(figname, bbox_inches='tight')

if __name__ == '__main__':
    # run("Psal.cds","aaa")
    rscu("Psal.RSCU.csv",'RSCU.pdf')