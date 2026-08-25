# %%
############### IMPORTING NEEDED MODULES FOR THE SCRIPT ###############
from operator import index
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import statistics
import sys
#from HNSC_immunescore import get_corr_pvals_IS
sys.path.append("E:/Immune_Evasion_Project/Mario/Scripts_Python")

from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_validate
from sklearn.metrics import recall_score
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import roc_auc_score, f1_score
import sklearn.metrics as metrics
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Ridge, RidgeCV, Lasso, LassoCV
from sklearn.preprocessing import scale
# %%
############### IMPORTING FUNCTIONS LIBRARY ###############

# We will import the 'Functions' library from our directory, make sure the file is present in the defined Directory above.
import Functions_Mario
from Functions_Mario import *

import Functions
from Functions import *

# This is for reloading the Library above. 
# import importlib
# importlib.reload(Functions_Mario)
# %%
############### IMPORTING REQUIRED DATA ###############

#HNSC_model = pd.read_csv("E:/Immune_Evasion_Project/Output/HNSC_model.csv", index_col=0)
HNSC_model = pd.read_csv("E:/Immune_Evasion_Project/Mario/Output/HNSC_model.csv", index_col=0)
HNSC_allparams =  pd.read_csv("E:/Immune_Evasion_Project/Mario/Cancer_Types/HNSC/HNSC_allgenes_parameters.csv", index_col=0)
Top100 = pd.read_csv("E:/Immune_Evasion_Project/Mario/Output/Top100_HNSC_metadata.csv", index_col=0)
pos_controls = Top100[Top100['Label'] == 1]
neg_controls = Top100[Top100['Label'] == 0]
# %%
# prepare data for logreg
neg_only = neg_controls
pos_only = pos_controls
mix_x = pd.concat([neg_only.iloc[:, 0:5],pos_only.iloc[:, 0:5]])
mix_y = pd.concat([neg_only['Label'],pos_only['Label']])

mix = pd.concat([neg_only, pos_only])
# %%
y = mix_y
x = mix_x #.iloc[:,[0,1,2,3,4]]
# %%
# FOR CROSS VAL
# getting train, test split for all data
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.2, random_state=42, stratify = y)
# %%

# fit the model with ALL variables and calculate auc score
clf_logistic = LogisticRegression(penalty = "none", max_iter=10000).fit(x_train,y_train)
preds = clf_logistic.predict_proba(x_test)[:,1] # only get prob of output being 1. 
pred = clf_logistic.predict(x_test)
#Using predict_proba vs predict will allow us to set a threshold to call it as 1 or 0

print(roc_auc_score(y_test, preds)) # choosing predict_proba instead of predict since we want to change the threshold to >0.5 to say its positive
# these scores can only be calculated using binary predictions and not pred probabilities
p = precision_score(y_test, pred)
r = recall_score(y_test, pred)
f = f1_score(y_test, pred)
print(f)
# %%

### USING ROC AUC ####
fpr, tpr, thresholds = metrics.roc_curve(y_test, preds)
auc = metrics.auc(fpr, tpr)
# %%
# calculate the g-mean for each threshold (sensitivity * specificity)
from numpy import sqrt, argmax,argmin
gmeans = sqrt(tpr * (1-fpr))
# locate the index of the largest g-mean
ix = argmax(gmeans)
print('Best Threshold=%f, G-Mean=%.3f' % (thresholds[ix], gmeans[ix]))
# plot the roc curve for the model
#plt.plot([0,1], [0,1], linestyle='--', label='No Skill')
sns.set(font_scale=1)
plt.plot(fpr, tpr, marker='.', label='Logistic')
plt.scatter(fpr[ix], tpr[ix], marker='o', color='black', label='Best')
# axis labels
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.legend()
# show the plot
plt.show()
# %%
plt.plot(fpr, tpr, label = 'AUC_predictproba= %0.7f' % auc)
plt.title('Logistic regression ROC')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive rate')
plt.legend(loc='best')
# %%
#### USING PRECISION RECALL CURVE ####
from sklearn.metrics import (precision_recall_curve, PrecisionRecallDisplay)

precision, recall, thresholds = precision_recall_curve(y_test, preds)
# plot the roc curve for the model
no_skill = len(y_test[y_test==1]) / len(y_test)
plt.plot([0,1], [no_skill,no_skill], linestyle='--', label='No Skill')
plt.plot(recall, precision, marker='.', label='Logistic')
# axis labels
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.legend()
# show the plot
plt.show()
# %%
# finding the best threshold using precision and recall (using F-score)
# convert to f score
fscore = (2 * precision * recall) / (precision + recall)
# locate the index of the largest f score
ix = argmax(fscore)
print('Best Threshold=%f, F-Score=%.3f' % (thresholds[ix], fscore[ix]))
# plot the roc curve for the model
no_skill = len(y_test[y_test==1]) / len(y_test)
plt.plot([0,1], [no_skill,no_skill], linestyle='--', label='No Skill')
plt.plot(recall, precision, marker='.', label='Logistic')
plt.scatter(recall[ix], precision[ix], marker='o', color='black', label='Best')
# axis labels
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.legend()
# show the plot
plt.show()

# %%
# define thresholds
from numpy import arange
thresholds = arange(0, 1, 0.001)
# apply threshold to positive probabilities to create labels
def to_labels(pos_probs, threshold):
	return (pos_probs >= threshold).astype('int')
# evaluate each threshold
scores = [f1_score(y_test, to_labels(preds, t)) for t in thresholds]
# get best threshold
ix = argmax(scores)
print('Threshold=%.3f, F-Score=%.5f' % (thresholds[ix], scores[ix]))
# %%
import statsmodels.api as sm

#define response variable
yy= y_train
#define predictor variables
xx = x_train
#add constant to predictor variables
xx = sm.add_constant(xx)

#fit linear regression model
log_reg = sm.Logit(yy, xx).fit()
#view model summary
print(log_reg.summary())
# %%
# LASSO
from numpy import arange

alphas = 10**np.linspace(10,-2,100)*0.5 # arange(0, 0.01, 0.001) 

lasso = LogisticRegression(penalty="l1",
    solver="liblinear",)
#Lasso(max_iter = 10000, normalize = True)
coefs = []

for a in alphas:
    lasso.set_params(C=a) # C is Inverse of regularization strength; smaller values specify stronger regularization.
    lasso.fit(scale(x_train), y_train)
    coefs.append(lasso.coef_.ravel().copy())
# %%
fig = plt.figure(figsize=(10,7))
ax = plt.gca()
ax.plot(alphas, coefs)
ax.set_xscale('log')
plt.axis('tight')
plt.xlabel('C')
plt.ylabel('weights')
plt.legend([ "Correlation SCNA and IS", "Correlation RNA and IS", "Correlation SCNA and RNA" ,"Z score using SCNA", "Z score using RNA expression"], loc="best")
#["SCNA_IS_unfiltered", "freqdel", "SCNA_RNA", "RNA_IS", "CCLE_SCNA_RNA", "CCLE_SCNA_prot", "survival_loss", "survival_gain", "survival_exp"])
#[ "corr SCNA and IS", "freq of deletion", "corr SCNA and RNA", "corr RNA and IS", "survival z score_loss","survival z score_gain", "survival z scorel_exp"])

#roc_auc_score(y_test, preds)

# %%
lassocv = LogisticRegressionCV(Cs=10, cv=3,penalty="l1", scoring="f1", solver="liblinear")
lassocv.fit(x_train, y_train) # solver already chooses the best alpha

lassocv_roc = LogisticRegressionCV(Cs=10, cv=3,penalty="l1", scoring="roc_auc", solver="liblinear")
lassocv_roc.fit(x_train, y_train) 


print("The best C chosen is", lassocv.C_)
print("coefficients:", lassocv.coef_)
print("coefficients_roc:", lassocv_roc.coef_)

print("roc_auc_score:",roc_auc_score(lassocv_roc.predict(x_test), y_test))
print("f1 score:",f1_score(lassocv.predict(x_test), y_test))
# %%
# CORRELATION MATRIX
corr_matrix = mix_x.corr(method = "spearman").round(2) # get corr, then round up to 2 decimal
mask = np.zeros_like(corr_matrix) # build matrix - get an array of zeros with same shape as corr_matrix  
mask[np.triu_indices_from(mask)] = True # mask (hide) the upper triangle
fig, ax = plt.subplots(figsize=(35,20))
sns.set(font_scale=4)
cmap = sns.diverging_palette(150, 275, as_cmap=True)
sns.heatmap(corr_matrix, mask = mask, annot=True, # The Seaborn heatmap ‘mask’ argument comes in handy when we want to cover part of the heatmap.
            cmap= cmap, vmin=-1, vmax=1)

# %%
## PREDICTIONS ##
all_params_rmNA = HNSC_model.dropna()
# create df for inference
real_preds = clf_logistic.predict_proba(all_params_rmNA)[:,1]

# %%
#first create dict
preds_dict = {}
for i in range(0,len(real_preds)):
    dict_vals = {all_params_rmNA.index[i]: real_preds[i]}
    preds_dict.update(dict_vals)

# then df
realpreds_df = pd.DataFrame(index = all_params_rmNA.index, columns = ["predict_prob", "class_outcome"])
realpreds_df["predict_prob"] = preds_dict.values()

for i in range(0,len(realpreds_df)):
    if realpreds_df.predict_prob[i] > 0.9:
        realpreds_df.class_outcome[i] = 1
    else:
       realpreds_df.class_outcome[i] = 0
# %%
realpreds_pos = realpreds_df[realpreds_df.class_outcome == 1] # predicted positives
realpreds_params = HNSC_model[HNSC_model.index.isin(realpreds_df.index)]
realpredspos_params = HNSC_model[HNSC_model.index.isin(realpreds_pos.index)]
# %%
realpreds_pos
realpreds_pos.shape
# %%
realpreds_params
realpreds_params.shape
# %%
realpredspos_params
realpredspos_params.shape
# %%
