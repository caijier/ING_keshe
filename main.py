import numpy as np
import sys
from matplotlib import pyplot as plt
#  判断是否表示数字

from sklearn.datasets import make_biclusters
from sklearn.cluster import SpectralCoclustering
from sklearn.metrics import accuracy_score
from sklearn.metrics import consensus_score
import sklearn.metrics as sm
import sklearn.ensemble as se  # 集合算法模块
import sklearn.utils as su  # 打乱数据

np.set_printoptions(threshold=np.inf)


def is_number(num):
    try:
        float(num)
        return True
    except ValueError:
        return False


str_int = {}
stick = 0


def trans(inpu):
    global stick
    global str_int
    if inpu in str_int:
        return str_int[inpu]
    str_int[inpu] = stick
    stick += 1
    return str_int[inpu]


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
row_dict = {}  # 基因名到行索引的映射字典
column_dict = {}  # 条件名到列索引的映射字典
row_list = []  # 行索引到基因名的映射字典
column_list = []  # 列索引到基因名的映射列表
row_num = 0  # 全局行索引
column_num = 0
min_list = []
max_list = []
ave_list = []
matrix = np.zeros((1300, 60000), dtype=np.float16)
column_index = 0


#  process file1 file2 file4
for file in [file1, file2, file4]:
    tem_gene_list = []  # 文件内基因索引映射到全局列表
    first_line = file.readline().split()
    for gene in first_line[1:]:
        if gene not in row_dict:  # 没有出现过的基因名
            row_dict[gene] = row_num
            row_list.append(gene)
            row_num += 1
        tem_gene_list.append(row_dict[gene])
    lines = file.readlines()  # 读取剩下行
    for line in lines:
        word = line.split()
        sum_column = 0  # 该条件的数据和
        num_valid = 0  # 有效数据个数
        maxd = sys.float_info.min
        mind = sys.float_info.max
        for item in word[1:]:
            if item != "NA":
                if maxd < float(item):
                    maxd = float(item)
                if mind > float(item):
                    mind = float(item)
                sum_column += float(item)
                num_valid += 1
        if num_valid / (len(word) - 1) < 0.7:  # 有效数据比例判定
            continue
        ave = (sum_column / num_valid - mind) / (maxd - mind)  # 该条件的归一化平均值

        #  记录全局该条件的数据信息
        min_list.append(mind)
        max_list.append(maxd)
        ave_list.append(ave)

        column_list.append(word[0])
        column_dict[word[0]] = column_num
        row_index = 0  # 文件行指向索引
        for item in word[1:]:
            if item != "NA":
                matrix[tem_gene_list[row_index]][column_index] = 100*(float(item) - mind) / (maxd - mind)  # 归一化处理
         #   else:
          #      matrix[tem_gene_list[row_index]][column_index] = ave
            row_index += 1
        column_index += 1
    print(column_index)

#  end


#  process file7
tem_gene_list = []  # 文件内基因索引映射到全局列表
first_line = file7.readline().split()
for gene in first_line[1:]:
    if gene not in row_dict:
        row_dict[gene] = row_num
        row_list.append(gene)
        row_num += 1
    tem_gene_list.append(row_dict[gene])
lines = file7.readlines()
for line in lines:
    word = line.split()
    sum_column = 0  # 该条件的数据和
    num_valid = 0  # 有效数据个数
    maxd = sys.float_info.min
    mind = sys.float_info.max
    for item in word[1:]:
        if not is_number(item):
            continue
        if item != 0:
            if maxd < float(item):
                maxd = float(item)
            if mind > float(item):
                mind = float(item)
            sum_column += float(item)
            num_valid += 1
    if mind == maxd or num_valid / (len(word) - 1) < 0.7:  # 剔除不需要的条件（全都相同或者有效数据比例不够）
        continue
    ave = (sum_column / num_valid - mind) / (maxd - mind)
    min_list.append(mind)
    max_list.append(maxd)
    ave_list.append(ave)
    column_list.append(word[0])
    column_dict[word[0]] = column_num
    row_index = 0
    for item in word[1:]:
        if not is_number(item):
            continue
        if row_index >= len(tem_gene_list):
            break
        if item != "NA":
            matrix[tem_gene_list[row_index]][column_index] = 100*(float(item) - mind) / (maxd - mind)
     #   else:
     #       matrix[tem_gene_list[row_index]][column_index] = ave
        row_index += 1
    column_index += 1
#  file7 process end

# process file5
tem_gene_list = []  # 文件内基因索引映射到全局列表
first_line = file5.readline().split()
for gene in first_line[2:]:
    if gene not in row_dict:
        row_dict[gene] = row_num
        row_list.append(gene)
        row_num += 1
    tem_gene_list.append(row_dict[gene])
lines = file5.readlines()
for line in lines:
    word = line.split()
    sum_column = 0  # 该条件的数据和
    column_list.append(word[0])
    column_dict[word[0]] = column_num
    row_index = 0

    min_list.append(0)
    max_list.append(1)
    ave_list.append(0.5)
    for item in word[1:]:
        if item == -2:
            matrix[tem_gene_list[row_index]][column_index] = 1
        elif item == -1:
            matrix[tem_gene_list[row_index]][column_index] = 25
        elif item == 0:
            matrix[tem_gene_list[row_index]][column_index] = 50
        elif item == 1:
            matrix[tem_gene_list[row_index]][column_index] = 75
        elif item == 2:
            matrix[tem_gene_list[row_index]][column_index] = 100
        row_index += 1
    column_index += 1
print(column_index)


# process file8  没有归一化
unsta_min = column_index
first_line = file8.readline().split()
for condi in first_line:
    column_list.append(column_index)
    column_dict[condi] = column_index
    min_list.append(sys.float_info.max)
    max_list.append(sys.float_info.min)
    ave_list.append(0)
    column_index += 1
unsta_max = column_index - 1

lines = file8.readlines()
for line in lines:
    column_index = unsta_min
    word = line.split()
    gene = word[0]
    if gene not in row_dict:
        row_dict[gene] = row_num
        row_list.append(gene)
        row_num += 1
    row_index = row_dict[gene]  # 获取全局行号
    for item in word[1:]:
        value = 100*float(item)
        matrix[row_index][column_index] = value
        if value < min_list[column_index]:
            min_list[column_index] = value
        if value > max_list[column_index]:
            max_list[column_index] = value
        if value != 0:
            ave_list[column_index] += 1  # 此时ave_list复用表示非零元素个数
        column_index += 1
print(unsta_max)
print("row_num", row_num)

#  跑对角线双聚类 每个基因打上0.。4 的标签

model = SpectralCoclustering(n_clusters=10, random_state=0)
model.fit(matrix)
for i in range(len(row_dict)):
    print(i, '.', row_list[i], ':', model.row_labels_[i])
print(model.column_labels_)



fit_data = matrix[np.argsort(model.row_labels_)]
fit_data = fit_data[:, np.argsort(model.column_labels_)]


#  临床变量随机森林
con_num = file9.readline().split().__len__() - 1
print("con_num:", con_num)
lines = file9.readlines()
sam_num = lines.__len__()
print("sam_num:", sam_num)
x_train = np.empty(shape=(sam_num, con_num), dtype=np.int)
y_train = np.empty(sam_num, dtype=int)
row_index = 0
for line in lines:
    col_index = 0
    p = 1
    word = line.split()
    while col_index < 12:
        if col_index <= 9:
            print(word[p])
            x_train[row_index][col_index] = trans(word[p])
        elif col_index == 9:
            if word[p] == 'NA':
                x_train[row_index][col_index] = trans(word[p])
            else:
                x_train[row_index][col_index] = trans(word[p]+word[p+1])
                p += 1
        else:
            x_train[row_index][col_index] = trans(word[p])
        col_index += 1
        p += 1
    if not word[0] in row_dict:
        y_train[row_index] = 0
    else:
        y_train[row_index] = model.row_labels_[row_dict[word[0]]]
    row_index += 1
print(x_train)
print(y_train)


x, y = su.shuffle(x_train, y_train, random_state=7)  # 打乱样本
train_size = int(len(x) * 0.9)
train_x, test_x, train_y, test_y = x[:train_size], x[train_size:], y[:train_size], y[train_size:]
model = se.RandomForestClassifier(n_estimators=100)
model = model.fit(train_x, train_y)
# 模型测试
pred_test_y = model.predict(test_x)
fi = model.feature_importances_
print(pred_test_y)
print(test_y)
# 模型评估
print('bike_hour的r2_score得分：', accuracy_score(test_y, pred_test_y))
print(fi)
#  process end
print("--->")
print(len(row_dict))

out_file = open("data.txt", "w")
for i in range(row_num):
    for j in range(unsta_max):
        print(matrix[i][j], file=out_file, end=' ')
    print(file=out_file)
out_file.close()
print(len(str_int))
"""
d = set()
for k in file[:7]:
    s = k.readline()
    for item in s.split():
        d.add(item)
    print(column_num)
    column_num += 1
    print(s)
    k.close()
lines = file8.readlines()
for line in lines:
    d.add(line.split()[0])
lines = file9.readlines()
for line in lines:
    d.add(line.split()[0])
print(d.__len__())
"""


