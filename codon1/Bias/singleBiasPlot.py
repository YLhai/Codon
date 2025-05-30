import pandas as pd
import matplotlib.pyplot as plt
import numpy as np



def single_biasplot(input_4ATCG,figName):
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
