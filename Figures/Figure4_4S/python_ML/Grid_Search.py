# %%
############### IMPORTING NEEDED LIBRARIES AND MODULES FOR THE SCRIPT ###############
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
from sklearn.model_selection import learning_curve
from sklearn.metrics import roc_auc_score
# %%
############### DEFINING ALL FUNCTIONS USED IN THIS SCRIPT ###############

# define a function to calculate the AUC
def calculateAUC(trueLabels,predProbabilities):
    return roc_auc_score(trueLabels,predProbabilities)

# define a function to invert the values of the binomial labels
def invert_label(x):
    return np.abs(x - 1)
# %%
############### Definining Control List ###############
# HNSC_metadata = pd.read_csv("E:/Immune_Evasion_Project/Mario/Output/DataControls_Dubrot_50.csv")
HNSC_metadata = pd.read_csv("E:/Immune_Evasion_Project/Mario/Output/Danias_allmetadata.csv")
HNSC_metadata
# %%
############### DATA MANIPULATION: INVERTING LABEL ###############

# apply the function to the column of binomial labels
# HNSC_metadata['Label'] = HNSC_metadata['Label'].apply(invert_label)

# print the resulting data frame
print(HNSC_metadata)

# %%
############### DEFINING 'x' and 'y' SETS ###############
x = HNSC_metadata.iloc[:,0:6]
y = HNSC_metadata[['Gene.Symbol','Label']]

x = x.set_index('Gene.Symbol')
y = y.set_index('Gene.Symbol')
# %%
############### GRID SEARCH ###############

from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import accuracy_score, log_loss
from sklearn.model_selection import train_test_split
import numpy as np

# Split the data into training, validation, and test sets
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.20, random_state=42, stratify=y)
x_train, x_val, y_train, y_val = train_test_split(x_train, y_train, test_size=0.20, random_state=42, stratify=y_train)

# Scale the features using StandardScaler
scaler = StandardScaler()
x_train = scaler.fit_transform(x_train)
x_val = scaler.transform(x_val)
x_test = scaler.transform(x_test)

# Define the logistic regression model
lr = LogisticRegression()

# Define the range of hyperparameters to search over
param_grid = {'penalty': ['l1', 'l2'], 'C': [0.1, 1, 10], 'solver': ['newton-cg', 'lbfgs', 'liblinear', 'sag', 'saga']}

# Perform grid search with cross-validation
grid = GridSearchCV(lr, param_grid, cv=5)
grid.fit(x_train, y_train)

# Get the best hyperparameters and corresponding performance
best_params = grid.best_params_
best_score = grid.best_score_
print('Best parameters:', best_params)
print('Best validation accuracy:', best_score)

# Evaluate the final model on the test set using the best hyperparameters
lr_final = LogisticRegression(**best_params)
lr_final.fit(x_train, y_train)
y_test_pred_label = lr_final.predict(x_test) 
y_test_pred_prob = lr_final.predict_proba(x_test)[:,1]
test_accuracy = accuracy_score(y_test, y_test_pred_label)
test_loss = log_loss(y_test, y_test_pred_prob)
print('Test accuracy:', test_accuracy)
print('Test log loss:', test_loss)