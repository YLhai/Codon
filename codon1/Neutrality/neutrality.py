import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from .GC123 import GC123

def neutrality_plot(GCstat,figname):
    # 生成示例数据
    data = pd.read_csv(GCstat, sep='\t')
    df = pd.DataFrame(data)
    df['GC12'] = (df['GC1'] + df['GC2']) / 2

    # 计算回归线（注意交换变量顺序）
    slope, intercept, r_value, p_value, std_err = stats.linregress(df['GC3'], df['GC12'])
    line = slope * df['GC3'] + intercept

    # 绘制中性图
    plt.figure(figsize=(8, 6))
    plt.scatter(df['GC3'], df['GC12'], c='steelblue', edgecolor='w', s=80, label='Genes')
    plt.plot(df['GC3'], line, '--', color='firebrick',
             label=f'Fit: y={slope:.4f}x+{intercept:.4f}\nR²={r_value**2:.4f}')
    plt.plot([20, 100], [30, 80], 'k:', lw=1, label='Neutral Line (y=x)')

    # 设置坐标轴和标签
    plt.xlabel('GC3 (%)', fontsize=12)
    plt.ylabel('GC12 (%)', fontsize=12)
    plt.title('Neutrality Plot (GC3 vs GC12)', pad=15, fontsize=14)
    plt.legend()
    plt.grid(alpha=0.3)
    plt.xlim(20, 100)
    plt.ylim(30, 80)
    plt.tight_layout()
    # 保存PDF
    plt.savefig(figname,
                format='pdf',
                dpi=300,
                bbox_inches='tight')

def run(cds_input,figname,outputPath,figType = "pdf"):

    out_results = outputPath + "GC123.stat"
    GC123(cds_input, out_results)

    inputstat = out_results
    figname = figname + "." + figType
    neutrality_plot(inputstat, figname)
