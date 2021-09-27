#import stuff
import numpy as np
import pandas as pd
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

def alleleFreqMatrix(AD, DP, fillna = True):
    #takes sparse AD and DP, returns dense AF matrix for plotting
    AD_df = pd.DataFrame(AD.todense())
    DP_df = pd.DataFrame(DP.todense())
    AF_df = AD_df/DP_df

    if fillna:
        AF_df = AF_df.fillna(0)

    return AF_df

def findKnee(BIC, sens=3):
    #Wrapper function for knee point locator given a series of deltaBIC

    #Remove negative BICs first

    BIC = BIC[BIC > 0]

    #Remove outliers first (Q3 + 1.5 IQR)
    q1 = np.percentile(BIC, 25) 
    q3 = np.percentile(BIC, 75)
    iqr = q3-q1
    t = q3 + 1.5*iqr ##threshold to determine outliers
    #print(t)
    ## if t is too small (ie. < 10k), set t to 10k
    if t < 10000: t = 10000

    ## remove outliers
    filtered_BIC = BIC[BIC <= t]

    y = np.sort(filtered_BIC.astype(float))
    x = np.linspace(0, 1, len(filtered_BIC)+1)[1:]
    kl = KneeLocator(x, y, curve="convex", direction="increasing", S=sens)

    #print(kl.knee_y)

    return x, y, kl.knee, kl.knee_y

if __name__ == '__main__':
    test = pd.read_csv('test/BIC_params.csv')
    findKnee(test.deltaBIC)