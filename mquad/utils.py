import os
from os import path
import sys
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from kneed import KneeLocator

def confusionMatrix(predicted_clone, real_label):

    conf_df = pd.DataFrame(data={'vireo': predicted_clone, 'real_label': real_label})
    confusion_matrix = pd.crosstab(conf_df['vireo'], conf_df['real_label'], rownames=['Predicted'], colnames=['Actual'])

    #of those cases predicted to belong to class c, which fraction truly belongs to class c? 
    precision = np.mean(confusion_matrix.max(axis=1)/confusion_matrix.sum(axis=1)) 

    #proportion of cases correctly identified as belonging to class c among all cases that truly belong to class c
    recall = np.mean(confusion_matrix.max(axis=0)/confusion_matrix.sum(axis=0)) 

    print('Precision = ' + str(precision))
    print('Recall = ' + str(recall))

    return confusion_matrix

def plot_confusionMatrix(mat, ax, cmap = 'Blues'):

    width, height = np.array(mat).shape
    text_colors = ['black', 'white']
    norm_conf = []

    for i in np.array(mat):
        a = 0
        tmp_arr = []
        a = sum(i, 0)
        for j in i:
            tmp_arr.append(float(j)/float(a))
        norm_conf.append(tmp_arr)

    res = ax.imshow(np.array(norm_conf), cmap=cmap, 
                    interpolation='nearest')
    
    for x in range(width):
        for y in range(height):
            ax.annotate(str(np.array(mat)[x][y]), xy=(y, x), 
                        horizontalalignment='center',
                        verticalalignment='center', color=text_colors[int(norm_conf[x][y] > 0.5)])

    return res