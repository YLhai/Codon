import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def glob_biasplot(input_ATCG,figName):

    # 读取数据
    # df = pd.read_csv('ATCG3', sep="\t")
    df = pd.read_csv(input_ATCG, sep="\t")
    # 分离平均值数据
    avg_row = df[df["AminoAcid"] == "Average"]
    df = df[df["AminoAcid"] != "Average"]

    # 定义可视化参数
    color_palette = ['red', 'orange', 'gold', 'green', 'cyan', 'blue', 'purple', 'black']  # 自定义颜色序列
    marker_size = 80  # 点大小


    # 计算坐标值
    def calc_ratio(g, c):
        return g / (g + c) if (g + c) > 0 else np.nan


    df["G_ratio"] = df.apply(lambda x: calc_ratio(x["G%"], x["C%"]), axis=1)
    df["A_ratio"] = df.apply(lambda x: calc_ratio(x["A%"], x["T%"]), axis=1)

    # 计算平均值坐标
    avg_g_ratio = calc_ratio(avg_row["G%"].values[0], avg_row["C%"].values[0])
    avg_a_ratio = calc_ratio(avg_row["A%"].values[0], avg_row["T%"].values[0])

    # 创建画布
    plt.figure(figsize=(10, 10))

    # 绘制每个氨基酸点（统一圆形+颜色序列）
    legend_handles = []
    for idx, (_, row) in enumerate(df.iterrows()):
        # 循环使用颜色序列
        color = color_palette[idx % len(color_palette)]

        scatter = plt.scatter(
            row["G_ratio"],
            row["A_ratio"],
            c=[color],
            marker='o',  # 统一使用圆形
            s=marker_size,
            edgecolor="w",
            linewidth=1,
            label=row["AminoAcid"]
        )
        legend_handles.append(scatter)

    # 添加平均值标记（红色×）
    avg_scatter = plt.scatter(
        avg_g_ratio,
        avg_a_ratio,
        marker="x",
        s=100,
        color="red",
        linewidths=2,
        label="Average"
    )
    legend_handles.append(avg_scatter)

    # 添加中性线
    plt.axhline(0.5, color="gray", linestyle="--", linewidth=1, alpha=0.7)
    plt.axvline(0.5, color="gray", linestyle="--", linewidth=1, alpha=0.7)

    # 设置坐标轴
    plt.xlabel("G3/(G3+C3)|4", fontsize=12)
    plt.ylabel("A3/(A3+T3)|4", fontsize=12)
    plt.title("Amino Acid Codon Usage Bias Plot", fontsize=14, pad=15)
    plt.xlim(0, 1)
    plt.ylim(0, 1)

    # 创建自定义图例
    plt.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor=(1, 1),
        frameon=True,
        framealpha=0.9,
        title="Amino Acids",
        handletextpad=0.5  # 调整图例项间距
    )

    # 调整布局
    plt.tight_layout()

    # 保存输出
    # plt.savefig("5.4Pcit.PR2-Bias.pdf",
    #             format='pdf',
    #             dpi=300,
    #             bbox_inches='tight')
    plt.savefig(figName,
                format='pdf',
                dpi=300,
                bbox_inches='tight')


def single_biasplot(input_4ATCG, figName):
    # 生成示例数据
    data = pd.read_csv(input_4ATCG)
    df = pd.DataFrame(data)

    # 计算PR2-Bias值
    df['A3_T3_ratio'] = df['A%'] / (df['A%'] + df['T%']).replace(0, np.nan)  # 避免除零错误
    df['G3_C3_ratio'] = df['G%'] / (df['G%'] + df['C%']).replace(0, np.nan)
    df = df.dropna()  # 移除无效值

    # 创建画布
    plt.figure(figsize=(10, 10))

    # 绘制主散点图
    main_scatter = plt.scatter(
        df['G3_C3_ratio'],
        df['A3_T3_ratio'],
        c='steelblue',
        alpha=0.6,
        edgecolor='w',
        s=80
    )

    # 添加密度等高线 (修正参数)
    # sns.kdeplot(
    #     data=df,  # 必须指定data参数
    #     x='A3_T3_ratio',
    #     y='G3_C3_ratio',
    #     levels=5,
    #     color='darkred',
    #     alpha=0.5,
    #     label='Density Contour'
    # )

    # 添加中性线
    plt.axhline(0.5, color='gray', linestyle='--', linewidth=1, alpha=0.7)
    plt.axvline(0.5, color='gray', linestyle='--', linewidth=1, alpha=0.7)

    # 设置坐标轴
    plt.xlabel('G3/(G3+C3)', fontsize=14)
    plt.ylabel('A3/(A3+T3)', fontsize=14)
    plt.title('PR2-Bias Plot Analysis with Density Contours', pad=20, fontsize=16)
    plt.xlim(0, 1)
    plt.ylim(0, 1)

    # 添加图例和网格
    # plt.legend()
    # plt.grid(linestyle='--', alpha=0.4)

    # 保存输出
    plt.tight_layout()
    plt.savefig(figName,
                format='pdf',
                dpi=300,
                bbox_inches='tight')

