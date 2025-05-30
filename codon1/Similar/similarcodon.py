import pandas as pd
import numpy as np
from ..Optimalcodon import RSCUofCDS



def filter_RSCU(RSCUfile):
    # 定义需要排除的密码子列表
    excluded_codons = ['ATG', 'TGG', 'TAA', 'TAG', 'TGA']
    df = pd.read_csv(RSCUfile, sep='\t')
    df_filter = df[~df['Codon'].isin(excluded_codons)]
    return df_filter

def calculate_similarity(df1, df2):

    merged = pd.merge(df1, df2, on='Codon', suffixes=('_1', '_2'))
    # 提取数值向量
    a = merged['RSCU_1'].values
    b = merged['RSCU_2'].values

    # 计算相似性
    if len(a) == 0 or len(b) == 0:
        return 0.0

    dot_product = np.dot(a, b)
    norm_a = np.linalg.norm(a)
    norm_b = np.linalg.norm(b)

    if norm_a == 0 or norm_b == 0:
        return 0.0

    R = dot_product / (norm_a * norm_b)
    D = (1-R)/2
    return D

def run(cdsA, cdsB,output_file):
    A_RSCU = output_file + 'A_RSCU.csv'
    B_RSCU = output_file + 'B_RSCU.csv'
    RSCUofCDS.run(cdsA, A_RSCU)
    RSCUofCDS.run(cdsB, B_RSCU)

    A = filter_RSCU(A_RSCU)
    B = filter_RSCU(B_RSCU)

    similarVlaue = calculate_similarity(A, B)
    print("SimilarValue:", similarVlaue)


if __name__ == "__main__":

    # Pcit = pd.read_csv('Pcit.RSCU.csv')
    # Post = pd.read_csv('Post.RSCU.csv')
    # Psal = pd.read_csv('Psal.RSCU.csv')
    # Pcit = filter_RSCU(Pcit)
    # Post = filter_RSCU(Post)
    # Psal = filter_RSCU(Psal)
    # Pcit_Post = calculate_similarity(Pcit, Post)
    # Pcit_Psal = calculate_similarity(Pcit, Psal)
    # Post_Psal = calculate_similarity(Post, Psal)
    #
    # print('Pcit_Post:')
    # print(Pcit_Post)
    # print('Pcit_Psal:')
    # print(Pcit_Psal)
    # print2
    pass
    # print(Post_Psal)
