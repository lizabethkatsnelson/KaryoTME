# %%
############### IMPORTING NEEDED MODULES FOR THE SCRIPT ###############
import pandas as pd
from sklearn.preprocessing import OneHotEncoder, LabelEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, precision_score, recall_score, f1_score, accuracy_score, confusion_matrix
# %%
############### IMPORTING REQUIRED DATA ###############

# Load the positive and negative control data
control_data = pd.read_csv('E:/Immune_Evasion_Project/Dania/Data/HNSC_Danias_Controls.csv')

# Load the 16k gene data
gene_data = pd.read_csv('E:/Immune_Evasion_Project/Dania/Data/HNSC_allgenes_parameters.csv')
# %%
############### DATA MANIPULATION ###############
# Remove the genes present in the control data from the gene data
gene_data = gene_data[~gene_data['Gene-Symbol'].isin(control_data['Gene-Symbol'])]

# Add a "population" label to the gene data
gene_data['class'] = "population"

# Merge the control and gene data
merged_data = pd.concat([control_data, gene_data]).reset_index(drop=True)

# Convert any numerical values in the class column to strings
merged_data['class'] = merged_data['class'].astype(str)
# %%
############### ENCODING ###############

# Encode the categorical labels as numerical values
label_encoder = LabelEncoder()
merged_data['class'] = label_encoder.fit_transform(merged_data['class'])
# %%
############### SPLITTING DATA ###############

# Split the data into training and test sets
X_train, X_test, y_train, y_test = train_test_split(merged_data_encoded.drop(['class'], axis=1), merged_data_encoded['class'], test_size=0.2, random_state=42)

# %%
# One-hot encode the categorical features
encoder = OneHotEncoder()
X_train = encoder.fit_transform(train_data[['Correlation SCNA and IS', 'Correlation RNA and IS', 'Correlation SCNA and RNA', 'Z score using SCNA', 'Z score using RNA expression']])


# Add the labels to the encoded data
merged_data_encoded['class'] = merged_data['class']
# %%
############### RUNNING LOGISTIC REGRESSION ###############
# Create the logistic regression model
model = LogisticRegression(multi_class='multinomial', solver='lbfgs')

# Fit the model to the training data
model.fit(X_train, y_train)

#%%
############### MODEL METRICS ###############
# Print the weights of the logistic regression model
print("Weights:", model.coef_)

# Predict the class labels for the test set
y_prob = model.predict_proba(X_test)

# Calculate the AUC score
auc = roc_auc_score(y_test, y_prob, multi_class='ovr')
print("AUC:", auc)

# Predict the class labels for the test set
y_pred = model.predict(X_test)

# Calculate the Precision score
precision = precision_score(y_test, y_pred, average='weighted')
print("Precision:", precision)

# Calculate the Recall score
recall = recall_score(y_test, y_pred, average='weighted')
print("Recall:", recall)

# Calculate the F1 score
f1 = f1_score(y_test, y_pred, average='weighted')
print("F1 Score:", f1)

# Calculate the overall accuracy of the model
accuracy = accuracy_score(y_test, y_pred)
print("Accuracy:", accuracy)

# Calculate the confusion matrix
conf_matrix = confusion_matrix(y_test, y_pred)
print("Confusion Matrix:\n", conf_matrix)
# %%
# Import libraries
import pandas as pd
from sklearn.preprocessing import OneHotEncoder, LabelEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, precision_score, recall_score, f1_score, accuracy_score, confusion_matrix

# Load the positive and negative control data
control_data = pd.read_csv('E:/Immune_Evasion_Project/Dania/Data/HNSC_Danias_Controls.csv')

# Load the 16k gene data
gene_data = pd.read_csv('E:/Immune_Evasion_Project/Dania/Data/HNSC_allgenes_parameters.csv')

# Remove the genes present in the control data from the gene data
gene_data = gene_data[~gene_data['Gene-Symbol'].isin(control_data['Gene-Symbol'])]

# Add a "population" label to the gene data
gene_data['class'] = "population"

# Merge the control and gene data
merged_data = pd.concat([control_data, gene_data]).reset_index(drop=True)

# Convert any numerical values in the class column to strings
merged_data['class'] = merged_data['class'].astype(str)

# Encode the categorical labels as numerical values
label_encoder = LabelEncoder()
merged_data['class'] = label_encoder.fit_transform(merged_data['class'])

# Split the data into training and test sets
train_data, test_data = train_test_split(merged_data, test_size=0.2, random_state=42)

# One-hot encode the categorical features in the training data
encoder = OneHotEncoder()
X_train = encoder.fit_transform(train_data[['Correlation SCNA and IS', 'Correlation RNA and IS', 'Correlation SCNA and RNA', 'Z score using SCNA', 'Z score using RNA expression']])

# Add the labels to the encoded training data
y_train = train_data['class']

# One-hot encode the categorical features in the test data
X_test = encoder.transform(test_data[['Correlation SCNA and IS', 'Correlation RNA and IS', 'Correlation SCNA and RNA', 'Z score using SCNA', 'Z score using RNA expression']])

# Convert the labels in the test set to categorical values
y_test = label_encoder.transform(test_data['class'])

# Create the logistic regression model
model = LogisticRegression(multi_class='multinomial', solver='lbfgs')

# Fit the model to the training data
model.fit(X_train, y_train)

# Print the weights of the logistic regression model
print("Weights:", model.coef_)

# Predict the class probabilities for the test set
y_prob = model.predict_proba(X_test)

# Calculate the AUC score
auc = roc_auc_score(y_test, y_prob, multi_class='ovr')
print("AUC:", auc)

# Predict the class labels for the test set
y_pred = model.predict(X_test)

# Calculate the Precision score
precision = precision_score(y_test, y_pred, average='weighted')
print("Precision:", precision)

# Calculate the Recall score
recall = recall_score(y_test, y_pred, average='weighted')
print("Recall:", recall)

# Calculate the F1 score
f1 = f1_score(y_test, y_pred, average='weighted')
print("F1 Score:", f1)

# Calculate the overall accuracy of the model
accuracy = accuracy_score(y_test, y_pred)
print("Accuracy:", accuracy)

# Calculate the confusion matrix
conf_matrix = confusion_matrix(y_test, y_pred)
print("Confusion Matrix:\n", conf_matrix)

# %%
