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
# Random Forest
from sklearn.ensemble import RandomForestClassifier
# %%
############### DEFINING 'x' and 'y' SETS ###############
controls = pd.read_csv("E:/Immune_Evasion_Project/Mario/Controls/MariosControls_models/MarioControls_Young_KIRC_MC38_100.csv")

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

# Initialize the logistic regression model without regularization
lr = LogisticRegression(penalty='none', max_iter=1000, random_state=123)

# Train the model on the training set
lr.fit(x_train, y_train)

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

# Perform 5-fold cross-validation on the Logistic Regression model
cv_scores = cross_val_score(lr, x_train, y_train, cv=5, scoring='accuracy')

# Compute the mean cross-validation accuracy
mean_cv_accuracy = np.mean(cv_scores)

# Compute the learning curves
train_sizes, train_scores, test_scores = learning_curve(lr, x_train, y_train, cv=5, train_sizes=np.linspace(0.1, 1.0, 10), scoring='accuracy')

# Plot the learning curves
plt.figure()
plt.title("Learning Curves (Logistic Regression)")
plt.xlabel("Training Examples")
plt.ylabel("Accuracy")
plt.plot(train_sizes, np.mean(train_scores, axis=1), 'o-', color="r", label="Training Accuracy")
plt.plot(train_sizes, np.mean(test_scores, axis=1), 'o-', color="g", label="Test Accuracy")
plt.legend(loc="best")
plt.show()

# Print the results
print("Training Accuracy:", train_accuracy)
print("Test Accuracy:", test_accuracy)
print("Training Loss:", train_loss)
print("Test Loss:", test_loss)

# Print the mean cross-validation accuracy
print("Mean Cross-Validation Accuracy (Logistic Regression):", mean_cv_accuracy)

# %%
############### LASSO REGRESSION ###############

# Initialize the logistic regression model with L1 regularization
lr_lasso = LogisticRegression(penalty='l1', solver='liblinear', C=1.0, max_iter=1000, random_state=123)

# Train the model on the training set
lr_lasso.fit(x_train, y_train)

# Predict the labels and probabilities on the training and test sets
y_train_pred_label = lr_lasso.predict(x_train)
y_train_pred_prob = lr_lasso.predict_proba(x_train)[:, 1]
y_test_pred_label = lr_lasso.predict(x_test)
y_test_pred_prob = lr_lasso.predict_proba(x_test)[:, 1]

# Compute the training and test accuracy
train_accuracy = accuracy_score(y_train, y_train_pred_label)
test_accuracy = accuracy_score(y_test, y_test_pred_label)

# Compute the training and test loss
train_loss = log_loss(y_train, y_train_pred_prob)
test_loss = log_loss(y_test, y_test_pred_prob)

# Perform 5-fold cross-validation on the Lasso Regression model
cv_scores = cross_val_score(lr_lasso, x_train, y_train, cv=5, scoring='accuracy')

# Compute the mean cross-validation accuracy
mean_cv_accuracy = np.mean(cv_scores)

# Compute the learning curves
train_sizes, train_scores, test_scores = learning_curve(lr_lasso, x_train, y_train, cv=5, train_sizes=np.linspace(0.1, 1.0, 10), scoring='accuracy')

# Plot the learning curves
plt.figure()
plt.title("Learning Curves (Lasso Regression)")
plt.xlabel("Training Examples")
plt.ylabel("Accuracy")
plt.plot(train_sizes, np.mean(train_scores, axis=1), 'o-', color="r", label="Training Accuracy")
plt.plot(train_sizes, np.mean(test_scores, axis=1), 'o-', color="g", label="Test Accuracy")
plt.legend(loc="best")
plt.show()

# Print the results
print("Training Accuracy:", train_accuracy)
print("Test Accuracy:", test_accuracy)
print("Training Loss:", train_loss)
print("Test Loss:", test_loss)

# Print the mean cross-validation accuracy
print("Mean Cross-Validation Accuracy (Lasso Regression):", mean_cv_accuracy)
# %%
############### RIDGE REGRESSION ###############

# Initialize the logistic regression model with L2 regularization
lr_ridge = LogisticRegression(penalty='l2', solver='lbfgs', C=1.0, max_iter=1000, random_state=123)

# Train the model on the training set
lr_ridge.fit(x_train, y_train)

# Predict the labels and probabilities on the training and test sets
y_train_pred_label = lr_ridge.predict(x_train)
y_train_pred_prob = lr_ridge.predict_proba(x_train)[:, 1]
y_test_pred_label = lr_ridge.predict(x_test)
y_test_pred_prob = lr_ridge.predict_proba(x_test)[:, 1]

# Compute the training and test accuracy
train_accuracy = accuracy_score(y_train, y_train_pred_label)
test_accuracy = accuracy_score(y_test, y_test_pred_label)

# Compute the training and test loss
train_loss = log_loss(y_train, y_train_pred_prob)
test_loss = log_loss(y_test, y_test_pred_prob)

# Perform 5-fold cross-validation on the Ridge Regression model
cv_scores = cross_val_score(lr_ridge, x_train, y_train, cv=5, scoring='accuracy')

# Compute the mean cross-validation accuracy
mean_cv_accuracy = np.mean(cv_scores)

# Compute the learning curves
train_sizes, train_scores, test_scores = learning_curve(lr_ridge, x_train, y_train, cv=5, train_sizes=np.linspace(0.1, 1.0, 10), scoring='accuracy')

# Plot the learning curves
plt.figure()
plt.title("Learning Curves (Ridge Regression)")
plt.xlabel("Training Examples")
plt.ylabel("Accuracy")
plt.plot(train_sizes, np.mean(train_scores, axis=1), 'o-', color="r", label="Training Accuracy")
plt.plot(train_sizes, np.mean(test_scores, axis=1), 'o-', color="g", label="Test Accuracy")
plt.legend(loc="best")
plt.show()

# Print the results
print("Training Accuracy:", train_accuracy)
print("Test Accuracy:", test_accuracy)
print("Training Loss:", train_loss)
print("Test Loss:", test_loss)

# Print the mean cross-validation accuracy
print("Mean Cross-Validation Accuracy (Ridge Regression):", mean_cv_accuracy)

# %%
############### K-NEAREST NEIGHBORS ###############

# Initialize the k-NN model with a fixed value for n_neighbors (e.g., 5)
knn = KNeighborsClassifier(n_neighbors=5)

# Train the model on the training set
knn.fit(x_train, y_train)

# Perform 5-fold cross-validation on the KNN model
cv_scores_knn = cross_val_score(knn, x_train, y_train, cv=5, scoring='accuracy')

# Compute the mean cross-validation accuracy
mean_cv_accuracy_knn = np.mean(cv_scores_knn)

# Predict the labels and probabilities on the training and test sets
y_train_pred_knn = knn.predict(x_train)
y_train_pred_knn_prob = knn.predict_proba(x_train)[:, 1]
y_test_pred_knn = knn.predict(x_test)
y_test_pred_knn_prob = knn.predict_proba(x_test)[:, 1]

# Compute the training and test accuracy
train_accuracy_knn = accuracy_score(y_train, y_train_pred_knn)
test_accuracy_knn = accuracy_score(y_test, y_test_pred_knn)

# Compute the training and test loss
train_loss_knn = log_loss(y_train, y_train_pred_knn_prob)
test_loss_knn = log_loss(y_test, y_test_pred_knn_prob)

# Print the results
print("Training Accuracy:", train_accuracy_knn)
print("Test Accuracy:", test_accuracy_knn)
print("Training Loss:", train_loss_knn)
print("Test Loss:", test_loss_knn)
print("Mean Cross-Validation Accuracy (KNN):", mean_cv_accuracy_knn)
# %%
############### SUPPORT VECTOR MACHINES ###############

# Initialize the SVM model with a fixed value for C (e.g., 1.0)
svm = SVC(kernel='linear', C=1.0, probability=True)

# Train the model on the training set
svm.fit(x_train, y_train)

# Perform 5-fold cross-validation on the SVM model
cv_scores_svm = cross_val_score(svm, x_train, y_train, cv=5, scoring='accuracy')

# Compute the mean cross-validation accuracy
mean_cv_accuracy_svm = np.mean(cv_scores_svm)

# Predict the labels and probabilities on the training and test sets
y_train_pred_svm = svm.predict(x_train)
y_train_pred_svm_prob = svm.predict_proba(x_train)[:, 1]
y_test_pred_svm = svm.predict(x_test)
y_test_pred_svm_prob = svm.predict_proba(x_test)[:, 1]

# Compute the training and test accuracy
train_accuracy_svm = accuracy_score(y_train, y_train_pred_svm)
test_accuracy_svm = accuracy_score(y_test, y_test_pred_svm)

# Compute the training and test loss
train_loss_svm = log_loss(y_train, y_train_pred_svm_prob)
test_loss_svm = log_loss(y_test, y_test_pred_svm_prob)

# Print the results
print("Training Accuracy:", train_accuracy_svm)
print("Test Accuracy:", test_accuracy_svm)
print("Training Loss:", train_loss_svm)
print("Test Loss:", test_loss_svm)
print("Mean Cross-Validation Accuracy (SVM):", mean_cv_accuracy_svm)
# %%
############### DECISION TREES ###############

# Initialize the decision tree model with a fixed value for max_depth (e.g., None)
dt = DecisionTreeClassifier(max_depth=None, random_state=123)

# Train the model on the training set
dt.fit(x_train, y_train)

# Perform 5-fold cross-validation on the decision tree model
cv_scores_dt = cross_val_score(dt, x_train, y_train, cv=5, scoring='accuracy')

# Compute the mean cross-validation accuracy
mean_cv_accuracy_dt = np.mean(cv_scores_dt)

# Predict the labels and probabilities on the training and test sets
y_train_pred_dt = dt.predict(x_train)
y_train_pred_dt_prob = dt.predict_proba(x_train)[:, 1]
y_test_pred_dt = dt.predict(x_test)
y_test_pred_dt_prob = dt.predict_proba(x_test)[:, 1]

# Compute the training and test accuracy
train_accuracy_dt = accuracy_score(y_train, y_train_pred_dt)
test_accuracy_dt = accuracy_score(y_test, y_test_pred_dt)

# Compute the training and test loss
train_loss_dt = log_loss(y_train, y_train_pred_dt_prob)
test_loss_dt = log_loss(y_test, y_test_pred_dt_prob)

# Print the results
print("Training Accuracy:", train_accuracy_dt)
print("Test Accuracy:", test_accuracy_dt)
print("Training Loss:", train_loss_dt)
print("Test Loss:", test_loss_dt)
print("Mean Cross-Validation Accuracy (Decision Tree):", mean_cv_accuracy_dt)
# %%
############### RANDOM FOREST ###############
# Initialize the Random Forest classifier
rf = RandomForestClassifier(random_state=123)

# Train the model on the training set
rf.fit(x_train, y_train)

# Predict the labels and probabilities on the training and test sets
y_train_pred_rf_label = rf.predict(x_train) 
y_train_pred_rf_prob = rf.predict_proba(x_train)[:,1] 
y_test_pred_rf_label = rf.predict(x_test) 
y_test_pred_rf_prob = rf.predict_proba(x_test)[:,1] 

# Compute the training and test accuracy
train_accuracy_rf = accuracy_score(y_train, y_train_pred_rf_label)
test_accuracy_rf = accuracy_score(y_test, y_test_pred_rf_label)

# Compute the training and test loss
train_loss_rf = log_loss(y_train, y_train_pred_rf_prob)
test_loss_rf = log_loss(y_test, y_test_pred_rf_prob)

# Perform 5-fold cross-validation on the Random Forest model
cv_scores_rf = cross_val_score(rf, x_train, y_train, cv=5, scoring='accuracy')

# Compute the mean cross-validation accuracy
mean_cv_accuracy_rf = np.mean(cv_scores_rf)

# Compute the learning curves
train_sizes_rf, train_scores_rf, test_scores_rf = learning_curve(rf, x_train, y_train, cv=5, train_sizes=np.linspace(0.1, 1.0, 10), scoring='accuracy')

# Plot the learning curves
plt.figure()
plt.title("Learning Curves (Random Forest)")
plt.xlabel("Training Examples")
plt.ylabel("Accuracy")
plt.plot(train_sizes_rf, np.mean(train_scores_rf, axis=1), 'o-', color="r", label="Training Accuracy")
plt.plot(train_sizes_rf, np.mean(test_scores_rf, axis=1), 'o-', color="g", label="Test Accuracy")
plt.legend(loc="best")
plt.show()

# Print the results
print("Training Accuracy:", train_accuracy_rf)
print("Test Accuracy:", test_accuracy_rf)
print("Training Loss:", train_loss_rf)
print("Test Loss:", test_loss_rf)

# Print the mean cross-validation accuracy
print("Mean Cross-Validation Accuracy (Random Forest):", mean_cv_accuracy_rf)
# %%