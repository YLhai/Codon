import sys

import pandas as pd
import matplotlib.pyplot as plt
import prince
from .RSCUofSingleGene import process_cds_file

def load_rscu_data(csv_path):
    df = pd.read_csv(csv_path, index_col=0)
    df = df.loc[:, (df != 0).any(axis=0)]
    return df


def perform_ca(df, n_components=59):
    ca = prince.CA(n_components=n_components, engine='sklearn')
    return ca.fit(df)


def plot_simple_ca(ca, df, output_path):
    # 获取前两个维度坐标
    row_coords = ca.row_coordinates(df).iloc[:, :2]

    # 创建画布
    fig, ax = plt.subplots(figsize=(10, 8))

    # 绘制基因点（不显示标签）
    ax.scatter(
        row_coords.iloc[:, 0],
        row_coords.iloc[:, 1],
        c='steelblue',
        alpha=0.6,
        s=30,
        edgecolor='w',
        linewidth=0.5
    )

    # 增强坐标轴显示
    ax.axhline(0, color='black', lw=0.8, ls='--')
    ax.axvline(0, color='black', lw=0.8, ls='--')
    ax.spines['bottom'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_linewidth(1.2)
    ax.spines['left'].set_linewidth(1.2)

    # 移除上右边框
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)

    # 设置坐标轴标签
    exp_var = ca.explained_inertia_
    ax.set_xlabel(f"CA1 ({exp_var[0] * 100:.1f}%)", fontsize=12)
    ax.set_ylabel(f"CA2 ({exp_var[1] * 100:.1f}%)", fontsize=12)

    # 移除网格和其他装饰
    ax.grid(False)
    ax.set_facecolor('white')

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def print_explained_variance(ca):
    print("前20个维度的解释比例：")
    print("------------------------")
    for i, var in enumerate(ca.explained_inertia_[:20]):
        print(f"维度 {i + 1:>2}: {var * 100:6.2f}%")
    print("------------------------")
    print(f"累计解释方差: {sum(ca.explained_inertia_[:20]) * 100:.2f}%")

def run(input_cds,figname,output_path,figType="pdf"):

    output_csv = output_path + "cas.rscu.csv"
    process_cds_file(input_cds, output_csv)
    input_csv = output_csv
    # 执行分析
    rscu_df = load_rscu_data(input_csv)
    ca_model = perform_ca(rscu_df)
    # 生成可视化
    figname = figname + "." + figType
    plot_simple_ca(ca_model, rscu_df, figname)

    # 打印解释比例
    print_explained_variance(ca_model)
    print(f"\n可视化结果已保存至：{figname}")


if __name__ == "__main__":
    # 参数设置



    input_cds = sys.argv[1]  # 输入文件路径
    figname = sys.argv[2]
    output_csv = "cas.rscu.csv"
    process_cds_file(input_cds, output_csv)

    input_csv = output_csv
    # 执行分析
    rscu_df = load_rscu_data(input_csv)
    ca_model = perform_ca(rscu_df)

    # 生成可视化
    plot_simple_ca(ca_model, rscu_df, figname)

    # 打印解释比例
    print_explained_variance(ca_model)
    print(f"\n可视化结果已保存至：{figname}")