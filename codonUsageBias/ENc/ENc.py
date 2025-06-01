import pandas as pd
import matplotlib.pyplot as plt
import numpy as np


def enc_GC3s(codonWout,figname,figType = "pdf"):
    # 示例数据（包含一个ENC>70的测试点）
    data = pd.read_csv(codonWout, sep='\t')
    df_row = pd.DataFrame(data)

    # 转换 Nc 列为数值，无效值设为 NaN
    df_row['Nc'] = pd.to_numeric(df_row['Nc'], errors='coerce')

    # 删除包含 NaN 的行
    df = df_row.dropna(subset=['Nc'])


    # 计算理论ENC曲线（注意避免除零错误）
    gc3s_range = np.linspace(0.01, 0.99, 100)
    enc_expected = 2 + gc3s_range + 29 / (gc3s_range**2 + (1 - gc3s_range)**2)

    # 创建画布
    plt.figure(figsize=(10, 8))

    # 绘制理论曲线（红色虚线）
    plt.plot(gc3s_range*100, enc_expected, 'r--', lw=2,
             label='Theoretical ENC (Mutation Bias Only)')

    # 绘制过滤后的数据点（蓝色散点）
    plt.scatter(df['GC3s']*100, df['Nc'],
                c='steelblue', edgecolor='w', s=80
                )

    # 设置坐标轴范围和标签
    plt.xlabel('GC3s (%)', fontsize=12)
    plt.ylabel('Effective Number of Codons (ENC)', fontsize=12)
    plt.title('ENC-GC3s Plot', fontsize=14, pad=20)
    plt.xlim(0, 100)    # X轴保持合理GC3s范围
    plt.ylim(0, 80)     # Y轴强制设置为0-80


    # 添加图例和网格
    plt.legend(loc='upper right')
    plt.grid(linestyle='--', alpha=0.4)

    # 显示和保存
    plt.tight_layout()
    plt.savefig(figname,
                format=figType,
                dpi=300,
                bbox_inches='tight')

def run(codonWout,figname,figType = "pdf"):
    figname = figname + "." + figType
    enc_GC3s(codonWout,figname,figType)
