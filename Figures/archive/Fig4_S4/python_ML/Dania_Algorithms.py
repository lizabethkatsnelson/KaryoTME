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
from sklearn.model_selection import train_test_split, cross_val_score, learning_curve, cross_validate
from sklearn.metrics import recall_score,accuracy_score,log_loss,precision_score,roc_auc_score, f1_score

import sklearn.metrics as metrics
from sklearn.preprocessing import StandardScaler,scale, OneHotEncoder
from sklearn.linear_model import Ridge, RidgeCV, Lasso, LassoCV
from sklearn import preprocessing
from sklearn.model_selection import learning_curve

# K-Nearest Neighbors
from sklearn.neighbors import KNeighborsClassifier
# Support Vector Machines
from sklearn.svm import SVC
from sklearn.model_selection import GridSearchCV #Needed also for Decision Trees 
# Decision Trees
from sklearn.tree import DecisionTreeClassifier
# %%
############### DEFINING 'x' and 'y' SETS ###############
x = controls[['Corr_SCNA_IS','Corr_RNA_IS','Corr_SCNA_RNA','Z.score_Loss','Z.score_Exp']]
y = controls[['Label']]
############### DATA MANIPULATION: SPLITTING DATA ###############
# getting train, test split for all data
x_train, x_test, y_train, y_test = train_test_split(x, y, test_size=0.20,random_state=42, stratify=y)

# Scale the features using StandardScaler
scaler = StandardScaler()
x_train = scaler.fit_transform(x_train)
x_test = scaler.transform(x_test)
# %%
############### LOGISTIC REGRESSION ###############

# Initialize the logistic regression model with no regularization
lr = LogisticRegression(penalty='none', max_iter=1000)

# Train the model on the training set
lr.fit(x_train, y_train)

# Perform 5-fold cross-validation
cv_scores = cross_val_score(lr, x_train, y_train, cv=5)
mean_cv_accuracy = np.mean(cv_scores)

# Predict the labels and probabilities on the training and test sets
y_train_pred_label = lr.predict(x_train)
y_train_pred_prob = lr.predict_proba(x_train)[:, 1]
y_test_pred_label = lr.predict(x_test)
y_test_pred_prob = lr.predict_proba(x_test)[:, 1]

# Compute the training and test accuracy
train_accuracy = accuracy_score(y_train, y_train_pred_label)
test_accuracy = accuracy_score(y_test, y_test_pred_label)

# Compute the training and test loss
train_loss = log_loss(y_train, y_train_pred_prob)
test_loss = log_loss(y_test, y_test_pred_prob)

# Compute the learning curves
train_sizes, train_scores, test_scores = learning_curve(lr, x_train, y_train, cv=5, train_sizes=np.linspace(0.1, 1.0, 10), scoring='accuracy')

# Plot the learning curves
plt.figure()
plt.title("Learning Curves (Logistic Regression)")
plt.xlabel("Training Examples")
plt.ylabel("Accuracy")
plt.plot(train_sizes, np.mean(train_scores, axis=1), 'o-', color="r", label="Training Accuracy")
plt.plot(train_sizes, np.mean(test_scores, axis=1), 'o-', color="g", label="Cross-Validation Accuracy")
plt.legend(loc="best")
plt.show()
# %%
############### LASSO REGRESSION ###############

# Initialize the logistic regression model with L1 regularization
lrcv = LogisticRegressionCV(penalty='l1', solver='liblinear', cv=5, max_iter=1000, random_state=123)

# Train the model on the training set
lrcv.fit(x_train, y_train)

# Predict the labels and probabilities on the training and test sets
y_train_pred_l1_label = lrcv.predict(x_train) 
y_train_pred_l1_prob = lrcv.predict_proba(x_train)[:,1] 
y_test_pred_l1_label = lrcv.predict(x_test) 
y_test_pred_l1_prob = lrcv.predict_proba(x_test)[:,1] 

# Compute the training and test accuracy
train_accuracy_l1 = accuracy_score(y_train, y_train_pred_l1_label)
test_accuracy_l1 = accuracy_score(y_test, y_test_pred_l1_label)

# Compute the training and test loss
train_loss_l1 = log_loss(y_train, y_train_pred_l1_prob)
test_loss_l1 = log_loss(y_test, y_test_pred_l1_prob)

# Compute the learning curves
train_sizes_l1, train_scores_l1, test_scores_l1 = learning_curve(lrcv, x_train, y_train, cv=5, train_sizes=np.linspace(0.1, 1.0, 10), scoring='accuracy')

# Plot the learning curves
plt.figure()
plt.title("Learning Curves (Lasso Regression)")
plt.xlabel("Training Examples")
plt.ylabel("Accuracy")
plt.plot(train_sizes_l1, np.mean(train_scores_l1, axis=1), 'o-', color="r", label="Training Accuracy")
plt.plot(train_sizes_l1, np.mean(test_scores_l1, axis=1), 'o-', color="g", label="Test Accuracy")
plt.legend(loc="best")
plt.show()

# Print the results
print("Training Accuracy:", train_accuracy_l1)
print("Test Accuracy:", test_accuracy_l1)
print("Training Loss:", train_loss_l1)
print("Test Loss:", test_loss_l1)
# %%
############### RIDGE REGRESSION ###############

# Define the range of C values to search over
Cs = [0.001, 0.01, 0.1, 1, 10]

# Define the logistic regression model with cross-validation
lr_cv = LogisticRegressionCV(Cs=Cs, cv=5, penalty='l2', max_iter=1000, solver='newton-cg')

# Fit the model on the training set
lr_cv.fit(x_train, y_train)

# Print the best hyperparameters
print('Best hyperparameters:', lr_cv.C_)

# Predict the labels and probabilities on the training and validation sets
y_train_pred_l2_label = lr_cv.predict(x_train) 
y_train_pred_l2_prob = lr_cv.predict_proba(x_train)[:,1] 

# Compute the training accuracy
train_accuracy_l2 = accuracy_score(y_train, y_train_pred_l2_label)

# Compute the training loss
train_loss_l2 = log_loss(y_train, y_train_pred_l2_prob)

# Predict the labels and probabilities on the test set using the model you just fitted above. 
y_test_pred_l2_label = lr_cv.predict(x_test) 
y_test_pred_l2_prob = lr_cv.predict_proba(x_test)[:,1]
# %% 
############### K-NEAREST NEIGHBORS ###############

# Define the range of K values to search over
K_values = list(range(1, 31))

# Initialize the k-NN model
knn = KNeighborsClassifier()

# Define the hyperparameter search space
param_grid = {'n_neighbors': K_values}

# Initialize the GridSearchCV object
grid_search = GridSearchCV(knn, param_grid, cv=5)

# Fit the model on the training set
grid_search.fit(x_train, y_train)

# Get the best K value
best_K = grid_search.best_params_['n_neighbors']

# Train the k-NN model with the best K value
best_knn = KNeighborsClassifier(n_neighbors=best_K)
best_knn.fit(x_train, y_train)

# Cross-validation
cv_scores_knn = cross_val_score(best_knn, x_train, y_train, cv=5)

# Predict the labels and probabilities on the test set
y_test_pred_knn = best_knn.predict(x_test)
y_test_pred_knn_prob = best_knn.predict_proba(x_test)[:, 1]
# %% 
############### SUPPORT VECTOR MACHINES ###############

from sklearn.svm import SVC
from sklearn.metrics import accuracy_score, roc_auc_score, precision_score, recall_score, f1_score
from sklearn.model_selection import GridSearchCV, cross_val_score

# Define the range of C values to search over
Cs = [0.001, 0.01, 0.1, 1, 10]

# Initialize the SVM model
svm = SVC(kernel='linear', probability=True)

# Define the hyperparameter search space
param_grid = {'C': Cs}

# Initialize the GridSearchCV object
grid_search = GridSearchCV(svm, param_grid, cv=5)

# Fit the model on the training set
grid_search.fit(x_train, y_train)

# Get the best C value
best_C = grid_search.best_params_['C']

# Train the SVM model with the best C value
best_svm = SVC(kernel='linear', C=best_C, probability=True)
best_svm.fit(x_train, y_train)

# Cross-validation
cv_scores_svm = cross_val_score(best_svm, x_train, y_train, cv=5)

# Predict the labels and probabilities on the test set
y_test_pred_svm = best_svm.predict(x_test)
y_test_pred_svm_prob = best_svm.predict_proba(x_test)[:, 1]
# %%
############### DECISION TREES ###############

# Define the range of max_depth values to search over
max_depths = [None, 3, 5, 7, 9]

# Initialize the Decision Tree model
dt = DecisionTreeClassifier(random_state=42)

# Create the grid search object
grid_dt = GridSearchCV(dt, param_grid={'max_depth': max_depths}, cv=5)

# Train the model on the training set
grid_dt.fit(x_train, y_train)

# Extract the best max_depth
best_max_depth_dt = grid_dt.best_params_['max_depth']

# Predict the labels and probabilities on the test set
y_test_pred_dt = grid_dt.predict(x_test)
y_test_pred_dt_prob = grid_dt.predict_proba(x_test)[:, 1]
# %%
############### LOGISTIC REGRESSION METRICS ###############
accuracies = []
aucs = []
precisions = []
recalls = []
f1_scores = []
weights = []

# Get the weights of the model
weights.append(lr.coef_)

# Get the accuracy, AUC, precision, recall, and f1 score on the test set
accuracy = accuracy_score(y_test, y_test_pred_label)
auc = roc_auc_score(y_test, y_test_pred_prob)
precision = precision_score(y_test, y_test_pred_label)
recall = recall_score(y_test, y_test_pred_label)
f1 = f1_score(y_test, y_test_pred_label)

# Append the results to the lists
accuracies.append(accuracy)
aucs.append(auc)
precisions.append(precision)
recalls.append(recall)
f1_scores.append(f1)

# Extracting Params
params1 = lr.get_params()

# Print the results
print('Max Iter Value:', params1['max_iter'])
print('Penalty Value:', params1['penalty'])
print('Random State Value:', params1['random_state'])
print('Solver Value:',params1['solver'])

print("Training Accuracy:", train_accuracy)
print("Training Loss:", train_loss)
print("Test Loss:", log_loss(y_test, y_test_pred_prob))
print("Test Accuracy:", np.mean(accuracies))
print("Test AUC:", np.mean(aucs))
print("Test Precision:", np.mean(precisions))
print("Test Recall:", np.mean(recalls))
print("Test F1 Score:", np.mean(f1_scores))

print("Cross-validated accuracy scores:", cv_scores)
print("Mean cross-validated accuracy:", mean_cv_accuracy)

print("Weights:", weights)
# %%
############### LASSO REGRESSION METRICS ###############

# Initialize the lists to store the results
accuracies_l1 = []
aucs_l1 = []
precisions_l1 = []
recalls_l1 = []
f1_scores_l1 = []
weights_l1 = []

# Get the weights of the model
weights_l1.append(lrcv.coef_)

# Get the predictions on the test data
y_test_pred_l1_label = lrcv.predict(x_test)
y_test_pred_l1_prob = lrcv.predict_proba(x_test)[:, 1]

# Get the accuracy, AUC, precision, recall, and f1 score on the test set
accuracy_l1 = accuracy_score(y_test, y_test_pred_l1_label)
auc_l1 = roc_auc_score(y_test, y_test_pred_l1_prob)
precision_l1 = precision_score(y_test, y_test_pred_l1_label)
recall_l1 = recall_score(y_test, y_test_pred_l1_label)
f1_l1 = f1_score(y_test, y_test_pred_l1_label)

# Append the results to the lists
accuracies_l1.append(accuracy_l1)
aucs_l1.append(auc_l1)
precisions_l1.append(precision_l1)
recalls_l1.append(recall_l1)
f1_scores_l1.append(f1_l1)

# Extracting Params
paramsl1 = lrcv.get_params()

# Print the results
print('CV Value:', paramsl1['cv'])
print('Max Iter Value:', paramsl1['max_iter'])
print('Penalty Value:', paramsl1['penalty'])
print('Random State Value:', paramsl1['random_state'])
print('Solver Value:',paramsl1['solver'])
print('Best hyperparameter C:', lrcv.C_)

print("Training Accuracy:", train_accuracy_l1)
print("Training Loss:", train_loss_l1)
print("Test Loss:", test_loss_l1)
print("Test Accuracy:", np.mean(accuracies_l1))
print("Test AUC:", np.mean(aucs_l1))
print("Test Precision:", np.mean(precisions_l1))
print("Test Recall:", np.mean(recalls_l1))
print("Test F1 Score:", np.mean(f1_scores_l1))

print("Weights:", weights_l1)
# %% 
############### RIDGE REGRESSION METRICS ###############

# Initialize the lists to store the results
accuracies_l2 = []
aucs_l2 = []
precisions_l2 = []
recalls_l2 = []
f1_scores_l2 = []
weights_l2 = []

# Get the weights of the model
weights_l2.append(lr_cv.coef_)
    
# Get the predictions on the test data
y_pred_l2 = lr_cv.predict(x_test)
    
# Get the accuracy, AUC, precision, recall, and f1 score on the test set
accuracy_l2 = accuracy_score(y_test, y_test_pred_l2_label)
auc_l2 = roc_auc_score(y_test, y_test_pred_l2_prob)
precision_l2 = precision_score(y_test, y_test_pred_l2_label)
recall_l2 = recall_score(y_test, y_test_pred_l2_label)
f1_l2 = f1_score(y_test, y_test_pred_l2_label)

# Append the results to the lists
accuracies_l2.append(accuracy_l2)
aucs_l2.append(auc_l2)
precisions_l2.append(precision_l2)
recalls_l2.append(recall_l2)
f1_scores_l2.append(f1_l2)

# Extracting Params
params = lr_cv.get_params()

# Print the results
print('CV Value:', params['cv'])
print('Max Iter Value:', params['max_iter'])
print('Penalty Value:', params['penalty'])
print('Random State Value:', params['random_state'])
print('Solver Value:', params['solver'])
print('Best hyperparameter C:', lr_cv.C_)

print("Training Accuracy:", train_accuracy_l2)
print("Test Accuracy:", np.mean(accuracies_l2))
print("Test AUC:", np.mean(aucs_l2))
print("Test Precision:", np.mean(precisions_l2))
print("Test Recall:", np.mean(recalls_l2))
print("Test F1 Score:", np.mean(f1_scores_l2))

print("Weights:", weights_l2)
# %%
############### K-NEAREST NEIGHBORS METRICS ###############

# Initialize the lists to store the results
accuracies_knn = []
aucs_knn = []
precisions_knn = []
recalls_knn = []
f1_scores_knn = []

# Get the accuracy, AUC, precision, recall, and f1 score on the test set
accuracy_knn = accuracy_score(y_test, y_test_pred_knn)
auc_knn = roc_auc_score(y_test, y_test_pred_knn_prob)
precision_knn = precision_score(y_test, y_test_pred_knn)
recall_knn = recall_score(y_test, y_test_pred_knn)
f1_knn = f1_score(y_test, y_test_pred_knn)

# Append the results to the lists
accuracies_knn.append(accuracy_knn)
aucs_knn.append(auc_knn)
precisions_knn.append(precision_knn)
recalls_knn.append(recall_knn)
f1_scores_knn.append(f1_knn)

# Print the results
print("Best K value:", best_K)
print("Test Accuracy:", np.mean(accuracies_knn))
print("Test AUC:", np.mean(aucs_knn))
print("Test Precision:", np.mean(precisions_knn))
print("Test Recall:", np.mean(recalls_knn))
print("Test F1 Score:", np.mean(f1_scores_knn))

print("Best K value:", best_K)
print("Cross-validated accuracy scores:", cv_scores_knn)
print("Mean cross-validated accuracy:", np.mean(cv_scores_knn))

# %%
############### SVM METRICS ###############

# Initialize the lists to store the results
accuracies_svm = []
aucs_svm = []
precisions_svm = []
recalls_svm = []
f1_scores_svm = []

# Get the accuracy, AUC, precision, recall, and f1 score on the test set
accuracy_svm = accuracy_score(y_test, y_test_pred_svm)
auc_svm = roc_auc_score(y_test, y_test_pred_svm_prob)
precision_svm = precision_score(y_test, y_test_pred_svm)
recall_svm = recall_score(y_test, y_test_pred_svm)
f1_svm = f1_score(y_test, y_test_pred_svm)

# Append the results to the lists
accuracies_svm.append(accuracy_svm)
aucs_svm.append(auc_svm)
precisions_svm.append(precision_svm)
recalls_svm.append(recall_svm)
f1_scores_svm.append(f1_svm)

# Print the results
print("Test Accuracy:", np.mean(accuracies_svm))
print("Test AUC:", np.mean(aucs_svm))
print("Test Precision:", np.mean(precisions_svm))
print("Test Recall:", np.mean(recalls_svm))
print("Test F1 Score:", np.mean(f1_scores_svm))

print("Best C value:", best_C)
print("Cross-validated accuracy scores:", cv_scores_svm)
print("Mean cross-validated accuracy:", np.mean(cv_scores_svm))

# %%
############### DECISION TREES METRICS ###############

# Initialize the lists to store the results
accuracies_dt = []
aucs_dt = []
precisions_dt = []
recalls_dt = []
f1_scores_dt = []
importances_dt = []

# Get the feature importances of the model
importances_dt.append(grid_dt.best_estimator_.feature_importances_)

# Get the accuracy, AUC, precision, recall, and f1 score on the test set
accuracy_dt = accuracy_score(y_test, y_test_pred_dt)
auc_dt = roc_auc_score(y_test, y_test_pred_dt_prob)
precision_dt = precision_score(y_test, y_test_pred_dt)
recall_dt = recall_score(y_test, y_test_pred_dt)
f1_dt = f1_score(y_test, y_test_pred_dt)

# Append the results to the lists
accuracies_dt.append(accuracy_dt)
aucs_dt.append(auc_dt)
precisions_dt.append(precision_dt)
recalls_dt.append(recall_dt)
f1_scores_dt.append(f1_dt)

# Print the results
print("Best max_depth:", best_max_depth_dt)
print("Test Accuracy:", np.mean(accuracies_dt))
print("Test AUC:", np.mean(aucs_dt))
print("Test Precision:", np.mean(precisions_dt))
print("Test Recall:", np.mean(recalls_dt))
print("Test F1 Score:", np.mean(f1_scores_dt))

print("Feature importances:", importances_dt)
# %%
############### CONFUSION MATRIX ###############

# Replace 'y_test' and 'y_test_pred_label' with the appropriate variable names for each algorithm
cm = confusion_matrix(y_test, y_test_pred_label)
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues')
plt.xlabel('Predicted')
plt.ylabel('True')
plt.title('Confusion Matrix')
plt.show()