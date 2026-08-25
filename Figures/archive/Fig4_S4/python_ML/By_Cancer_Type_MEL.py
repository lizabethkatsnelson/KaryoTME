#%%
### IMPORTING ALL NECESSARY MODULES

from operator import index
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
import statistics
import sys
sys.path.append("E:/Immune_Evasion_Project/Scripts")

# LogReg modules
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
from sklearn.preprocessing import OneHotEncoder
from sklearn import preprocessing

# %%
############################################################# BRCA_MEL #############################################################################

## Definining Control List
BRCA_MEL_metadata = pd.read_csv("E:/Immune_Evasion_Project/Output/DataControls_BRCA_AP2_MEL.csv")
BRCA_MEL_metadata
# %%
# Definining 'x' and 'y' sets
x = BRCA_MEL_metadata.iloc[:,0:6]
y = BRCA_MEL_metadata[['Gene_ID','Label']]

x = x.set_index('Gene_ID')
y = y.set_index('Gene_ID')
# y['Labels'] = y['Labels'].astype(np.float64)
# %%
# FOR CROSS VAL
# getting train, test split for all data
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.20,random_state=42, stratify=y)

# %%
############################################################# DEFINING ALGORITHMS #############################################################################

#################### LOGISTIC REGRESSION ####################

#We have provided a helper function to calculate the AUC using sklearn:
from sklearn.metrics import roc_auc_score
def calculateAUC(trueLabels,predProbabilities):
  return roc_auc_score(trueLabels,predProbabilities)
#We suggest using this function:
from sklearn.linear_model import LogisticRegression

#Fit a LogisticRegression model to the X_train data. 

#Use a random_state of 123 with no regularization/penalty.
lr = LogisticRegressionCV(cv=5, max_iter=10000,solver='liblinear').fit(x_train, y_train) 

#Predict the labels and probabilities on the test set using the model you just fitted above. 
y_test_pred_label = lr.predict(x_test) 
y_test_pred_prob = lr.predict_proba(x_test)[:,1] 

#################### LASSO REGRESSION ####################

#We have provided a helper function to calculate the AUC using sklearn:
from sklearn.metrics import roc_auc_score
def calculateAUC(trueLabels,predProbabilities):
  return roc_auc_score(trueLabels,predProbabilities)

#Fitting a LogisticRegressionCV model to the 'x_train' data. 
#Using a random_state of 123, with L1 penalty, 5-fold cross validation, and a Solver parameter
lrcv = LogisticRegressionCV(penalty = 'l1',cv=5,solver='liblinear',max_iter=10000).fit(x_train, y_train)

#Predicting the labels and probabilities on the test set using the model that was just fitted above. 
y_test_pred_l1_label = lrcv.predict(x_test) 
y_test_pred_l1_prob = lrcv.predict_proba(x_test)[:,1] 

# %%

#### BRCA — LOGISTIC REGRESSION METRICS ####
# Initialize the lists to store the results
accuracies = []
aucs = []
precisions = []
recalls = []
f1_scores = []
weights = []

# Get the weights of the model
weights.append(lr.coef_)

# Get the accuracy, AUC, precision, recall, and f1 score
accuracy = accuracy_score(y_test, y_test_pred_prob.round())
auc = calculateAUC(y_test, y_test_pred_prob)
precision = precision_score(y_test, y_test_pred_prob.round())
recall = recall_score(y_test, y_test_pred_prob.round())
f1 = f1_score(y_test, y_test_pred_prob.round())

# Append the results to the lists
accuracies.append(accuracy)
aucs.append(auc)
precisions.append(precision)
recalls.append(recall)
f1_scores.append(f1)
    
# Print the results
print("Weights:", weights)
print("Accuracy:", np.mean(accuracies))
print("AUC:", np.mean(aucs))
print("Precision:", np.mean(precisions))
print("Recall:", np.mean(recalls))
print("F1 Score:", np.mean(f1_scores))

# %%
#### BRCA — LASSO REGRESSION METRICS ####
# Initialize the lists to store the results
accuracies = []
aucs = []
precisions = []
recalls = []
f1_scores = []
weights = []

# Get the weights of the model
weights.append(lrcv.coef_)
    
# Get the predictions on the test data
y_pred = lrcv.predict(x_test)
    
# Get the accuracy, AUC, precision, recall, and f1 score
accuracy = accuracy_score(y_test, y_test_pred_l1_prob.round())
auc = calculateAUC(y_test, y_test_pred_l1_prob)
precision = precision_score(y_test, y_test_pred_l1_prob.round())
recall = recall_score(y_test, y_test_pred_l1_prob.round())
f1 = f1_score(y_test, y_test_pred_l1_prob.round())

# Append the results to the lists
accuracies.append(accuracy)
aucs.append(auc)
precisions.append(precision)
recalls.append(recall)
f1_scores.append(f1)
    
# Print the results
print("Weights:", weights)
print("Accuracy:", np.mean(accuracies))
print("AUC:", np.mean(aucs))
print("Precision:", np.mean(precisions))
print("Recall:", np.mean(recalls))
print("F1 Score:", np.mean(f1_scores))

# %%
############################################################# COADREAD_MEL #############################################################################

## Definining Control List
COADREAD_MEL_metadata = pd.read_csv("E:/Immune_Evasion_Project/Output/DataControls_COADREAD_AP2_MEL.csv")
COADREAD_MEL_metadata
# %%
#################### DEFINING 'x' and 'y' SETS ####################

x = COADREAD_MEL_metadata.iloc[:,0:6]
y = COADREAD_MEL_metadata[['Gene_ID','Label']]

x = x.set_index('Gene_ID')
y = y.set_index('Gene_ID')
# y['Labels'] = y['Labels'].astype(np.float64)
# %%
#################### SPLITTING ####################

# getting train, test split for all data
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.20,random_state=42, stratify=y)

# %%
############################################################# DEFINING ALGORITHMS #############################################################################

#################### LOGISTIC REGRESSION ####################

#We have provided a helper function to calculate the AUC using sklearn:
from sklearn.metrics import roc_auc_score
def calculateAUC(trueLabels,predProbabilities):
  return roc_auc_score(trueLabels,predProbabilities)
#We suggest using this function:
from sklearn.linear_model import LogisticRegression

#Fit a LogisticRegression model to the X_train data. 

#Use a random_state of 123 with no regularization/penalty.
lr = LogisticRegressionCV(cv=5, max_iter=10000,solver='liblinear').fit(x_train, y_train) 

#Predict the labels and probabilities on the test set using the model you just fitted above. 
y_test_pred_label = lr.predict(x_test) 
y_test_pred_prob = lr.predict_proba(x_test)[:,1] 

#################### LASSO REGRESSION ####################

#We have provided a helper function to calculate the AUC using sklearn:
from sklearn.metrics import roc_auc_score
def calculateAUC(trueLabels,predProbabilities):
  return roc_auc_score(trueLabels,predProbabilities)

#Fitting a LogisticRegressionCV model to the 'x_train' data. 
#Using a random_state of 123, with L1 penalty, 5-fold cross validation, and a Solver parameter
lrcv = LogisticRegressionCV(penalty = 'l1',cv=5,solver='liblinear',max_iter=10000).fit(x_train, y_train)

#Predicting the labels and probabilities on the test set using the model that was just fitted above. 
y_test_pred_l1_label = lrcv.predict(x_test) 
y_test_pred_l1_prob = lrcv.predict_proba(x_test)[:,1] 

# %%

#### COADREAD — LOGISTIC REGRESSION METRICS ####
# Initialize the lists to store the results
accuracies = []
aucs = []
precisions = []
recalls = []
f1_scores = []
weights = []

# Get the weights of the model
weights.append(lr.coef_)

# Get the accuracy, AUC, precision, recall, and f1 score
accuracy = accuracy_score(y_test, y_test_pred_prob.round())
auc = calculateAUC(y_test, y_test_pred_prob)
precision = precision_score(y_test, y_test_pred_prob.round())
recall = recall_score(y_test, y_test_pred_prob.round())
f1 = f1_score(y_test, y_test_pred_prob.round())

# Append the results to the lists
accuracies.append(accuracy)
aucs.append(auc)
precisions.append(precision)
recalls.append(recall)
f1_scores.append(f1)
    
# Print the results
print("Weights:", weights)
print("Accuracy:", np.mean(accuracies))
print("AUC:", np.mean(aucs))
print("Precision:", np.mean(precisions))
print("Recall:", np.mean(recalls))
print("F1 Score:", np.mean(f1_scores))

# %%
#### COADREAD — LASSO REGRESSION METRICS ####
# Initialize the lists to store the results
accuracies = []
aucs = []
precisions = []
recalls = []
f1_scores = []
weights = []

# Get the weights of the model
weights.append(lrcv.coef_)
    
# Get the predictions on the test data
y_pred = lrcv.predict(x_test)
    
# Get the accuracy, AUC, precision, recall, and f1 score
accuracy = accuracy_score(y_test, y_test_pred_l1_prob.round())
auc = calculateAUC(y_test, y_test_pred_l1_prob)
precision = precision_score(y_test, y_test_pred_l1_prob.round())
recall = recall_score(y_test, y_test_pred_l1_prob.round())
f1 = f1_score(y_test, y_test_pred_l1_prob.round())

# Append the results to the lists
accuracies.append(accuracy)
aucs.append(auc)
precisions.append(precision)
recalls.append(recall)
f1_scores.append(f1)
    
# Print the results
print("Weights:", weights)
print("Accuracy:", np.mean(accuracies))
print("AUC:", np.mean(aucs))
print("Precision:", np.mean(precisions))
print("Recall:", np.mean(recalls))
print("F1 Score:", np.mean(f1_scores))
# %%
############################################################# LUAD_MEL #############################################################################

## Definining Control List
LUAD_MEL_metadata = pd.read_csv("E:/Immune_Evasion_Project/Output/DataControls_LUAD_AP2_MEL.csv")
LUAD_MEL_metadata
# %%
#################### DEFINING 'x' and 'y' SETS ####################

x = LUAD_MEL_metadata.iloc[:,0:6]
y = LUAD_MEL_metadata[['Gene_ID','Label']]

x = x.set_index('Gene_ID')
y = y.set_index('Gene_ID')
# y['Labels'] = y['Labels'].astype(np.float64)
# %%
#################### SPLITTING ####################

# getting train, test split for all data
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.20,random_state=42, stratify=y)

# %%
############################################################# DEFINING ALGORITHMS #############################################################################

#################### LOGISTIC REGRESSION ####################

#We have provided a helper function to calculate the AUC using sklearn:
from sklearn.metrics import roc_auc_score
def calculateAUC(trueLabels,predProbabilities):
  return roc_auc_score(trueLabels,predProbabilities)
#We suggest using this function:
from sklearn.linear_model import LogisticRegression

#Fit a LogisticRegression model to the X_train data. 

#Use a random_state of 123 with no regularization/penalty.
lr = LogisticRegressionCV(cv=5, max_iter=10000,solver='liblinear').fit(x_train, y_train) 

#Predict the labels and probabilities on the test set using the model you just fitted above. 
y_test_pred_label = lr.predict(x_test) 
y_test_pred_prob = lr.predict_proba(x_test)[:,1] 

#################### LASSO REGRESSION ####################

#We have provided a helper function to calculate the AUC using sklearn:
from sklearn.metrics import roc_auc_score
def calculateAUC(trueLabels,predProbabilities):
  return roc_auc_score(trueLabels,predProbabilities)

#Fitting a LogisticRegressionCV model to the 'x_train' data. 
#Using a random_state of 123, with L1 penalty, 5-fold cross validation, and a Solver parameter
lrcv = LogisticRegressionCV(penalty = 'l1',cv=5,solver='liblinear',max_iter=10000).fit(x_train, y_train)

#Predicting the labels and probabilities on the test set using the model that was just fitted above. 
y_test_pred_l1_label = lrcv.predict(x_test) 
y_test_pred_l1_prob = lrcv.predict_proba(x_test)[:,1] 

# %%

#### LUAD — LOGISTIC REGRESSION METRICS ####
# Initialize the lists to store the results
accuracies = []
aucs = []
precisions = []
recalls = []
f1_scores = []
weights = []

# Get the weights of the model
weights.append(lr.coef_)

# Get the accuracy, AUC, precision, recall, and f1 score
accuracy = accuracy_score(y_test, y_test_pred_prob.round())
auc = calculateAUC(y_test, y_test_pred_prob)
precision = precision_score(y_test, y_test_pred_prob.round())
recall = recall_score(y_test, y_test_pred_prob.round())
f1 = f1_score(y_test, y_test_pred_prob.round())

# Append the results to the lists
accuracies.append(accuracy)
aucs.append(auc)
precisions.append(precision)
recalls.append(recall)
f1_scores.append(f1)
    
# Print the results
print("Weights:", weights)
print("Accuracy:", np.mean(accuracies))
print("AUC:", np.mean(aucs))
print("Precision:", np.mean(precisions))
print("Recall:", np.mean(recalls))
print("F1 Score:", np.mean(f1_scores))

# %%
#### LUAD — LASSO REGRESSION METRICS ####
# Initialize the lists to store the results
accuracies = []
aucs = []
precisions = []
recalls = []
f1_scores = []
weights = []

# Get the weights of the model
weights.append(lrcv.coef_)
    
# Get the predictions on the test data
y_pred = lrcv.predict(x_test)
    
# Get the accuracy, AUC, precision, recall, and f1 score
accuracy = accuracy_score(y_test, y_test_pred_l1_prob.round())
auc = calculateAUC(y_test, y_test_pred_l1_prob)
precision = precision_score(y_test, y_test_pred_l1_prob.round())
recall = recall_score(y_test, y_test_pred_l1_prob.round())
f1 = f1_score(y_test, y_test_pred_l1_prob.round())

# Append the results to the lists
accuracies.append(accuracy)
aucs.append(auc)
precisions.append(precision)
recalls.append(recall)
f1_scores.append(f1)
    
# Print the results
print("Weights:", weights)
print("Accuracy:", np.mean(accuracies))
print("AUC:", np.mean(aucs))
print("Precision:", np.mean(precisions))
print("Recall:", np.mean(recalls))
print("F1 Score:", np.mean(f1_scores))
# %%
############################################################# LUSC_MEL #############################################################################

## Definining Control List
LUSC_MEL_metadata = pd.read_csv("E:/Immune_Evasion_Project/Output/DataControls_LUSC_AP2_MEL.csv")
LUSC_MEL_metadata
# %%
#################### DEFINING 'x' and 'y' SETS ####################

x = LUSC_MEL_metadata.iloc[:,0:6]
y = LUSC_MEL_metadata[['Gene_ID','Label']]

x = x.set_index('Gene_ID')
y = y.set_index('Gene_ID')
# y['Labels'] = y['Labels'].astype(np.float64)
# %%
#################### SPLITTING ####################

# getting train, test split for all data
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.20,random_state=42, stratify=y)

# %%
############################################################# DEFINING ALGORITHMS #############################################################################

#################### LOGISTIC REGRESSION ####################

#We have provided a helper function to calculate the AUC using sklearn:
from sklearn.metrics import roc_auc_score
def calculateAUC(trueLabels,predProbabilities):
  return roc_auc_score(trueLabels,predProbabilities)
#We suggest using this function:
from sklearn.linear_model import LogisticRegression

#Fit a LogisticRegression model to the X_train data. 

#Use a random_state of 123 with no regularization/penalty.
lr = LogisticRegressionCV(cv=5, max_iter=10000,solver='liblinear').fit(x_train, y_train) 

#Predict the labels and probabilities on the test set using the model you just fitted above. 
y_test_pred_label = lr.predict(x_test) 
y_test_pred_prob = lr.predict_proba(x_test)[:,1] 

#################### LASSO REGRESSION ####################

#We have provided a helper function to calculate the AUC using sklearn:
from sklearn.metrics import roc_auc_score
def calculateAUC(trueLabels,predProbabilities):
  return roc_auc_score(trueLabels,predProbabilities)

#Fitting a LogisticRegressionCV model to the 'x_train' data. 
#Using a random_state of 123, with L1 penalty, 5-fold cross validation, and a Solver parameter
lrcv = LogisticRegressionCV(penalty = 'l1',cv=5,solver='liblinear',max_iter=10000).fit(x_train, y_train)

#Predicting the labels and probabilities on the test set using the model that was just fitted above. 
y_test_pred_l1_label = lrcv.predict(x_test) 
y_test_pred_l1_prob = lrcv.predict_proba(x_test)[:,1] 

# %%

#### LUSC — LOGISTIC REGRESSION METRICS ####
# Initialize the lists to store the results
accuracies = []
aucs = []
precisions = []
recalls = []
f1_scores = []
weights = []

# Get the weights of the model
weights.append(lr.coef_)

# Get the accuracy, AUC, precision, recall, and f1 score
accuracy = accuracy_score(y_test, y_test_pred_prob.round())
auc = calculateAUC(y_test, y_test_pred_prob)
precision = precision_score(y_test, y_test_pred_prob.round())
recall = recall_score(y_test, y_test_pred_prob.round())
f1 = f1_score(y_test, y_test_pred_prob.round())

# Append the results to the lists
accuracies.append(accuracy)
aucs.append(auc)
precisions.append(precision)
recalls.append(recall)
f1_scores.append(f1)
    
# Print the results
print("Weights:", weights)
print("Accuracy:", np.mean(accuracies))
print("AUC:", np.mean(aucs))
print("Precision:", np.mean(precisions))
print("Recall:", np.mean(recalls))
print("F1 Score:", np.mean(f1_scores))

# %%
#### LUSC — LASSO REGRESSION METRICS ####
# Initialize the lists to store the results
accuracies = []
aucs = []
precisions = []
recalls = []
f1_scores = []
weights = []

# Get the weights of the model
weights.append(lrcv.coef_)
    
# Get the predictions on the test data
y_pred = lrcv.predict(x_test)
    
# Get the accuracy, AUC, precision, recall, and f1 score
accuracy = accuracy_score(y_test, y_test_pred_l1_prob.round())
auc = calculateAUC(y_test, y_test_pred_l1_prob)
precision = precision_score(y_test, y_test_pred_l1_prob.round())
recall = recall_score(y_test, y_test_pred_l1_prob.round())
f1 = f1_score(y_test, y_test_pred_l1_prob.round())

# Append the results to the lists
accuracies.append(accuracy)
aucs.append(auc)
precisions.append(precision)
recalls.append(recall)
f1_scores.append(f1)
    
# Print the results
print("Weights:", weights)
print("Accuracy:", np.mean(accuracies))
print("AUC:", np.mean(aucs))
print("Precision:", np.mean(precisions))
print("Recall:", np.mean(recalls))
print("F1 Score:", np.mean(f1_scores))
# %%
############################################################# SKCM_MEL #############################################################################

## Definining Control List
SKCM_MEL_metadata = pd.read_csv("E:/Immune_Evasion_Project/Output/DataControls_SKCM_AP2_MEL.csv")
SKCM_MEL_metadata
# %%
#################### DEFINING 'x' and 'y' SETS ####################

x = SKCM_MEL_metadata.iloc[:,0:6]
y = SKCM_MEL_metadata[['Gene_ID','Label']]

x = x.set_index('Gene_ID')
y = y.set_index('Gene_ID')
# y['Labels'] = y['Labels'].astype(np.float64)
# %%
#################### SPLITTING ####################

# getting train, test split for all data
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.20,random_state=42, stratify=y)

# %%
############################################################# DEFINING ALGORITHMS #############################################################################

#################### LOGISTIC REGRESSION ####################

#We have provided a helper function to calculate the AUC using sklearn:
from sklearn.metrics import roc_auc_score
def calculateAUC(trueLabels,predProbabilities):
  return roc_auc_score(trueLabels,predProbabilities)
#We suggest using this function:
from sklearn.linear_model import LogisticRegression

#Fit a LogisticRegression model to the X_train data. 

#Use a random_state of 123 with no regularization/penalty.
lr = LogisticRegressionCV(cv=5, max_iter=10000,solver='liblinear').fit(x_train, y_train) 

#Predict the labels and probabilities on the test set using the model you just fitted above. 
y_test_pred_label = lr.predict(x_test) 
y_test_pred_prob = lr.predict_proba(x_test)[:,1] 

#################### LASSO REGRESSION ####################

#We have provided a helper function to calculate the AUC using sklearn:
from sklearn.metrics import roc_auc_score
def calculateAUC(trueLabels,predProbabilities):
  return roc_auc_score(trueLabels,predProbabilities)

#Fitting a LogisticRegressionCV model to the 'x_train' data. 
#Using a random_state of 123, with L1 penalty, 5-fold cross validation, and a Solver parameter
lrcv = LogisticRegressionCV(penalty = 'l1',cv=5,solver='liblinear',max_iter=10000).fit(x_train, y_train)

#Predicting the labels and probabilities on the test set using the model that was just fitted above. 
y_test_pred_l1_label = lrcv.predict(x_test) 
y_test_pred_l1_prob = lrcv.predict_proba(x_test)[:,1] 

# %%

#### SKCM — LOGISTIC REGRESSION METRICS ####
# Initialize the lists to store the results
accuracies = []
aucs = []
precisions = []
recalls = []
f1_scores = []
weights = []

# Get the weights of the model
weights.append(lr.coef_)

# Get the accuracy, AUC, precision, recall, and f1 score
accuracy = accuracy_score(y_test, y_test_pred_prob.round())
auc = calculateAUC(y_test, y_test_pred_prob)
precision = precision_score(y_test, y_test_pred_prob.round())
recall = recall_score(y_test, y_test_pred_prob.round())
f1 = f1_score(y_test, y_test_pred_prob.round())

# Append the results to the lists
accuracies.append(accuracy)
aucs.append(auc)
precisions.append(precision)
recalls.append(recall)
f1_scores.append(f1)
    
# Print the results
print("Weights:", weights)
print("Accuracy:", np.mean(accuracies))
print("AUC:", np.mean(aucs))
print("Precision:", np.mean(precisions))
print("Recall:", np.mean(recalls))
print("F1 Score:", np.mean(f1_scores))

# %%
#### SKCM — LASSO REGRESSION METRICS ####
# Initialize the lists to store the results
accuracies = []
aucs = []
precisions = []
recalls = []
f1_scores = []
weights = []

# Get the weights of the model
weights.append(lrcv.coef_)
    
# Get the predictions on the test data
y_pred = lrcv.predict(x_test)
    
# Get the accuracy, AUC, precision, recall, and f1 score
accuracy = accuracy_score(y_test, y_test_pred_l1_prob.round())
auc = calculateAUC(y_test, y_test_pred_l1_prob)
precision = precision_score(y_test, y_test_pred_l1_prob.round())
recall = recall_score(y_test, y_test_pred_l1_prob.round())
f1 = f1_score(y_test, y_test_pred_l1_prob.round())

# Append the results to the lists
accuracies.append(accuracy)
aucs.append(auc)
precisions.append(precision)
recalls.append(recall)
f1_scores.append(f1)
    
# Print the results
print("Weights:", weights)
print("Accuracy:", np.mean(accuracies))
print("AUC:", np.mean(aucs))
print("Precision:", np.mean(precisions))
print("Recall:", np.mean(recalls))
print("F1 Score:", np.mean(f1_scores))
# %%
############################################################# PAAD_MEL #############################################################################

## Definining Control List
PAAD_MEL_metadata = pd.read_csv("E:/Immune_Evasion_Project/Output/DataControls_PAAD_AP2_MEL.csv")
PAAD_MEL_metadata
# %%
#################### DEFINING 'x' and 'y' SETS ####################

x = PAAD_MEL_metadata.iloc[:,0:6]
y = PAAD_MEL_metadata[['Gene_ID','Label']]

x = x.set_index('Gene_ID')
y = y.set_index('Gene_ID')
# y['Labels'] = y['Labels'].astype(np.float64)
# %%
#################### SPLITTING ####################

# getting train, test split for all data
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.20,random_state=42, stratify=y)

# %%
############################################################# DEFINING ALGORITHMS #############################################################################

#################### LOGISTIC REGRESSION ####################

#We have provided a helper function to calculate the AUC using sklearn:
from sklearn.metrics import roc_auc_score
def calculateAUC(trueLabels,predProbabilities):
  return roc_auc_score(trueLabels,predProbabilities)
#We suggest using this function:
from sklearn.linear_model import LogisticRegression

#Fit a LogisticRegression model to the X_train data. 

#Use a random_state of 123 with no regularization/penalty.
lr = LogisticRegressionCV(cv=5, max_iter=10000,solver='liblinear').fit(x_train, y_train) 

#Predict the labels and probabilities on the test set using the model you just fitted above. 
y_test_pred_label = lr.predict(x_test) 
y_test_pred_prob = lr.predict_proba(x_test)[:,1] 

#################### LASSO REGRESSION ####################

#We have provided a helper function to calculate the AUC using sklearn:
from sklearn.metrics import roc_auc_score
def calculateAUC(trueLabels,predProbabilities):
  return roc_auc_score(trueLabels,predProbabilities)

#Fitting a LogisticRegressionCV model to the 'x_train' data. 
#Using a random_state of 123, with L1 penalty, 5-fold cross validation, and a Solver parameter
lrcv = LogisticRegressionCV(penalty = 'l1',cv=5,solver='liblinear',max_iter=10000).fit(x_train, y_train)

#Predicting the labels and probabilities on the test set using the model that was just fitted above. 
y_test_pred_l1_label = lrcv.predict(x_test) 
y_test_pred_l1_prob = lrcv.predict_proba(x_test)[:,1] 

# %%

#### PAAD — LOGISTIC REGRESSION METRICS ####
# Initialize the lists to store the results
accuracies = []
aucs = []
precisions = []
recalls = []
f1_scores = []
weights = []

# Get the weights of the model
weights.append(lr.coef_)

# Get the accuracy, AUC, precision, recall, and f1 score
accuracy = accuracy_score(y_test, y_test_pred_prob.round())
auc = calculateAUC(y_test, y_test_pred_prob)
precision = precision_score(y_test, y_test_pred_prob.round())
recall = recall_score(y_test, y_test_pred_prob.round())
f1 = f1_score(y_test, y_test_pred_prob.round())

# Append the results to the lists
accuracies.append(accuracy)
aucs.append(auc)
precisions.append(precision)
recalls.append(recall)
f1_scores.append(f1)
    
# Print the results
print("Weights:", weights)
print("Accuracy:", np.mean(accuracies))
print("AUC:", np.mean(aucs))
print("Precision:", np.mean(precisions))
print("Recall:", np.mean(recalls))
print("F1 Score:", np.mean(f1_scores))

# %%
#### PAAD — LASSO REGRESSION METRICS ####
# Initialize the lists to store the results
accuracies = []
aucs = []
precisions = []
recalls = []
f1_scores = []
weights = []

# Get the weights of the model
weights.append(lrcv.coef_)
    
# Get the predictions on the test data
y_pred = lrcv.predict(x_test)
    
# Get the accuracy, AUC, precision, recall, and f1 score
accuracy = accuracy_score(y_test, y_test_pred_l1_prob.round())
auc = calculateAUC(y_test, y_test_pred_l1_prob)
precision = precision_score(y_test, y_test_pred_l1_prob.round())
recall = recall_score(y_test, y_test_pred_l1_prob.round())
f1 = f1_score(y_test, y_test_pred_l1_prob.round())

# Append the results to the lists
accuracies.append(accuracy)
aucs.append(auc)
precisions.append(precision)
recalls.append(recall)
f1_scores.append(f1)
    
# Print the results
print("Weights:", weights)
print("Accuracy:", np.mean(accuracies))
print("AUC:", np.mean(aucs))
print("Precision:", np.mean(precisions))
print("Recall:", np.mean(recalls))
print("F1 Score:", np.mean(f1_scores))
# %%
############################################################# KIRC_MEL #############################################################################

## Definining Control List
KIRC_MELmetadata = pd.read_csv("E:/Immune_Evasion_Project/Output/DataControls_KIRC_AP2_MEL.csv")
KIRC_MEL_metadata
# %%
#################### DEFINING 'x' and 'y' SETS ####################

x = KIRC_MEL_metadata.iloc[:,0:6]
y = KIRC_MEL_metadata[['Gene_ID','Label']]

x = x.set_index('Gene_ID')
y = y.set_index('Gene_ID')
# y['Labels'] = y['Labels'].astype(np.float64)
# %%
#################### SPLITTING ####################

# getting train, test split for all data
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.20,random_state=42, stratify=y)

# %%
############################################################# DEFINING ALGORITHMS #############################################################################

#################### LOGISTIC REGRESSION ####################

#We have provided a helper function to calculate the AUC using sklearn:
from sklearn.metrics import roc_auc_score
def calculateAUC(trueLabels,predProbabilities):
  return roc_auc_score(trueLabels,predProbabilities)
#We suggest using this function:
from sklearn.linear_model import LogisticRegression

#Fit a LogisticRegression model to the X_train data. 

#Use a random_state of 123 with no regularization/penalty.
lr = LogisticRegressionCV(cv=5, max_iter=10000,solver='liblinear').fit(x_train, y_train) 

#Predict the labels and probabilities on the test set using the model you just fitted above. 
y_test_pred_label = lr.predict(x_test) 
y_test_pred_prob = lr.predict_proba(x_test)[:,1] 

#################### LASSO REGRESSION ####################

#We have provided a helper function to calculate the AUC using sklearn:
from sklearn.metrics import roc_auc_score
def calculateAUC(trueLabels,predProbabilities):
  return roc_auc_score(trueLabels,predProbabilities)

#Fitting a LogisticRegressionCV model to the 'x_train' data. 
#Using a random_state of 123, with L1 penalty, 5-fold cross validation, and a Solver parameter
lrcv = LogisticRegressionCV(penalty = 'l1',cv=5,solver='liblinear',max_iter=10000).fit(x_train, y_train)

#Predicting the labels and probabilities on the test set using the model that was just fitted above. 
y_test_pred_l1_label = lrcv.predict(x_test) 
y_test_pred_l1_prob = lrcv.predict_proba(x_test)[:,1] 

# %%

#### KIRC — LOGISTIC REGRESSION METRICS ####
# Initialize the lists to store the results
accuracies = []
aucs = []
precisions = []
recalls = []
f1_scores = []
weights = []

# Get the weights of the model
weights.append(lr.coef_)

# Get the accuracy, AUC, precision, recall, and f1 score
accuracy = accuracy_score(y_test, y_test_pred_prob.round())
auc = calculateAUC(y_test, y_test_pred_prob)
precision = precision_score(y_test, y_test_pred_prob.round())
recall = recall_score(y_test, y_test_pred_prob.round())
f1 = f1_score(y_test, y_test_pred_prob.round())

# Append the results to the lists
accuracies.append(accuracy)
aucs.append(auc)
precisions.append(precision)
recalls.append(recall)
f1_scores.append(f1)
    
# Print the results
print("Weights:", weights)
print("Accuracy:", np.mean(accuracies))
print("AUC:", np.mean(aucs))
print("Precision:", np.mean(precisions))
print("Recall:", np.mean(recalls))
print("F1 Score:", np.mean(f1_scores))

# %%
#### KIRC — LASSO REGRESSION METRICS ####
# Initialize the lists to store the results
accuracies = []
aucs = []
precisions = []
recalls = []
f1_scores = []
weights = []

# Get the weights of the model
weights.append(lrcv.coef_)
    
# Get the predictions on the test data
y_pred = lrcv.predict(x_test)
    
# Get the accuracy, AUC, precision, recall, and f1 score
accuracy = accuracy_score(y_test, y_test_pred_l1_prob.round())
auc = calculateAUC(y_test, y_test_pred_l1_prob)
precision = precision_score(y_test, y_test_pred_l1_prob.round())
recall = recall_score(y_test, y_test_pred_l1_prob.round())
f1 = f1_score(y_test, y_test_pred_l1_prob.round())

# Append the results to the lists
accuracies.append(accuracy)
aucs.append(auc)
precisions.append(precision)
recalls.append(recall)
f1_scores.append(f1)
    
# Print the results
print("Weights:", weights)
print("Accuracy:", np.mean(accuracies))
print("AUC:", np.mean(aucs))
print("Precision:", np.mean(precisions))
print("Recall:", np.mean(recalls))
print("F1 Score:", np.mean(f1_scores))
# %%
############################################################# KIRP_MEL #############################################################################

## Definining Control List
KIRP_MEL_metadata = pd.read_csv("E:/Immune_Evasion_Project/Output/DataControls_KIRP_AP2_MEL.csv")
KIRP_MEL_metadata
# %%
#################### DEFINING 'x' and 'y' SETS ####################

x = KIRP_MEL_metadata.iloc[:,0:6]
y = KIRP_MEL_metadata[['Gene_ID','Label']]

x = x.set_index('Gene_ID')
y = y.set_index('Gene_ID')
# y['Labels'] = y['Labels'].astype(np.float64)
# %%
#################### SPLITTING ####################

# getting train, test split for all data
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.20,random_state=42, stratify=y)

# %%
############################################################# DEFINING ALGORITHMS #############################################################################

#################### LOGISTIC REGRESSION ####################

#We have provided a helper function to calculate the AUC using sklearn:
from sklearn.metrics import roc_auc_score
def calculateAUC(trueLabels,predProbabilities):
  return roc_auc_score(trueLabels,predProbabilities)
#We suggest using this function:
from sklearn.linear_model import LogisticRegression

#Fit a LogisticRegression model to the X_train data. 

#Use a random_state of 123 with no regularization/penalty.
lr = LogisticRegressionCV(cv=5, max_iter=10000,solver='liblinear').fit(x_train, y_train) 

#Predict the labels and probabilities on the test set using the model you just fitted above. 
y_test_pred_label = lr.predict(x_test) 
y_test_pred_prob = lr.predict_proba(x_test)[:,1] 

#################### LASSO REGRESSION ####################

#We have provided a helper function to calculate the AUC using sklearn:
from sklearn.metrics import roc_auc_score
def calculateAUC(trueLabels,predProbabilities):
  return roc_auc_score(trueLabels,predProbabilities)

#Fitting a LogisticRegressionCV model to the 'x_train' data. 
#Using a random_state of 123, with L1 penalty, 5-fold cross validation, and a Solver parameter
lrcv = LogisticRegressionCV(penalty = 'l1',cv=5,solver='liblinear',max_iter=10000).fit(x_train, y_train)

#Predicting the labels and probabilities on the test set using the model that was just fitted above. 
y_test_pred_l1_label = lrcv.predict(x_test) 
y_test_pred_l1_prob = lrcv.predict_proba(x_test)[:,1] 

# %%

#### KIRP — LOGISTIC REGRESSION METRICS ####
# Initialize the lists to store the results
accuracies = []
aucs = []
precisions = []
recalls = []
f1_scores = []
weights = []

# Get the weights of the model
weights.append(lr.coef_)

# Get the accuracy, AUC, precision, recall, and f1 score
accuracy = accuracy_score(y_test, y_test_pred_prob.round())
auc = calculateAUC(y_test, y_test_pred_prob)
precision = precision_score(y_test, y_test_pred_prob.round())
recall = recall_score(y_test, y_test_pred_prob.round())
f1 = f1_score(y_test, y_test_pred_prob.round())

# Append the results to the lists
accuracies.append(accuracy)
aucs.append(auc)
precisions.append(precision)
recalls.append(recall)
f1_scores.append(f1)
    
# Print the results
print("Weights:", weights)
print("Accuracy:", np.mean(accuracies))
print("AUC:", np.mean(aucs))
print("Precision:", np.mean(precisions))
print("Recall:", np.mean(recalls))
print("F1 Score:", np.mean(f1_scores))

# %%
#### KIRP — LASSO REGRESSION METRICS ####
# Initialize the lists to store the results
accuracies = []
aucs = []
precisions = []
recalls = []
f1_scores = []
weights = []

# Get the weights of the model
weights.append(lrcv.coef_)
    
# Get the predictions on the test data
y_pred = lrcv.predict(x_test)
    
# Get the accuracy, AUC, precision, recall, and f1 score
accuracy = accuracy_score(y_test, y_test_pred_l1_prob.round())
auc = calculateAUC(y_test, y_test_pred_l1_prob)
precision = precision_score(y_test, y_test_pred_l1_prob.round())
recall = recall_score(y_test, y_test_pred_l1_prob.round())
f1 = f1_score(y_test, y_test_pred_l1_prob.round())

# Append the results to the lists
accuracies.append(accuracy)
aucs.append(auc)
precisions.append(precision)
recalls.append(recall)
f1_scores.append(f1)
    
# Print the results
print("Weights:", weights)
print("Accuracy:", np.mean(accuracies))
print("AUC:", np.mean(aucs))
print("Precision:", np.mean(precisions))
print("Recall:", np.mean(recalls))
print("F1 Score:", np.mean(f1_scores))
# %%
