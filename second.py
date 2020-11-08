import sklearn as sk
import numpy   as np

file1 = open(r'.\brca\miRNA'
             r'\miRNA_HiSeq_gene.txt', 'r', encoding='ANSI')  # mirna数据集

file2 = open(r'.\brca\蛋白质表达RPPA'
             r'\RPPA.txt', 'r', encoding='ANSI')  # 蛋白质表达

file3 = open(r'.\brca\甲基化HumanMethylation450'
             r'\HumanMethylation450', 'r', encoding='ANSI')  # 甲基化

file4 = open(r'.\brca\基因表达RNA_seq'
             r'\HiSeqV2_log.txt', 'r', encoding='ANSI')  # 基因表达

file5 = open(r'.\brca\拷贝数变异数据cnv'
             r'\Gistic2_CopyNumber_Gistic2_all_thresholded_by_genes.txt', 'r', encoding='ANSI')  # 拷贝数变异阈值型

file6 = open(r'.\brca\体细胞突变somatic mutation'
             r'\mutation_wustl_hiseq_gene.txt', 'r', encoding='ANSI')  # 体细胞突变

file7 = open(r'.\brca\通路活性Pathway activity Pathway_Paradigm_RNASeq'
             r'\Pathway_Paradigm_RNASeq.txt', 'r', encoding='ANSI')  # 通路活性

file8 = open(r'.\brca\转录因子调控影响Transcription factor regulatory impact - HiSeqV2, by RABIT'
             r'\RABIT_BRCA.HiSeq _HiSeqV2.txt', 'r', encoding='ANSI')  # 转录因子

file9 = open(r'.\brca\临床变量clinical'
             r'\BRCA_clinicalMatrix_1.txt', 'r', encoding='UTF-8')  # 临床变量

file = [file1, file2, file3, file4, file5, file6, file7, file8, file9]
num = 1
for k in file:
    print(num)
    num += 1
    print(k.readline())
    k.close()




