# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Replace the following dictionaries with the Test AUC values for each algorithm and cancer type
BRCA_auc = {'LogReg': 0.7448, 'Lasso': 0.7656, 'Ridge': 0.7447, 'KNN': 0, 'SVM': 0.7656, 'DT': 0.7447,'RandomForest': 0.7447}
COADREAD_auc = {'LogReg': 0.5260, 'Lasso': 0.5468, 'Ridge': 0.5260, 'KNN': 0, 'SVM': 0.5156, 'DT': 0.5260, 'RandomForest': 0.5260}
HNSC_auc = {'LogReg': 0.6145, 'Lasso': 0.6145, 'Ridge': 0.6145, 'KNN': 0, 'SVM': 0.6354, 'DT': 0.3072, 'RandomForest': 0.6145}
KIRC_auc = {'LogReg':0.6666, 'Lasso':  0.6875, 'Ridge': 0.6666, 'KNN': 0, 'SVM': 0.3645, 'DT': 0.3333, 'RandomForest': 0.6666}
KIRP_auc = {'LogReg': 0.6093, 'Lasso': 0.5937, 'Ridge': 0.6093, 'KNN': 0, 'SVM': 0.5937, 'DT': 0.3046, 'RandomForest': 0.6093}
LUAD_auc = {'LogReg': 0.6302, 'Lasso': 0.6406, 'Ridge': 0.6302, 'KNN': 0, 'SVM': 0.6406, 'DT': 0.3151, 'RandomForest': 0.6302}
LUSC_auc = {'LogReg': 0.5989, 'Lasso': 0.6197, 'Ridge': 0.5989, 'KNN': 0, 'SVM': 0.3854, 'DT': 0.2994, 'RandomForest': 0.5989}
PAAD_auc = {'LogReg': 0.5156, 'Lasso': 0.5052, 'Ridge': 0.5156, 'KNN': 0, 'SVM': 0.4687, 'DT': 0.2578, 'RandomForest': 0.5156}
SKCM_auc = {'LogReg': 0.6927, 'Lasso': 0.6979, 'Ridge': 0.6927, 'KNN': 0, 'SVM': 0.6614, 'DT': 0.3463, 'RandomForest': 0.6927}

# List of cancer types
cancer_types = ['BRCA', 'COADREAD', 'HNSC', 'KIRC', 'KIRP', 'LUAD', 'LUSC', 'PAAD', 'SKCM']

# List of Test AUC values dictionaries for each cancer type
auc_values = [BRCA_auc, COADREAD_auc, HNSC_auc, KIRC_auc, KIRP_auc, LUAD_auc, LUSC_auc, PAAD_auc, SKCM_auc]

# Create a DataFrame
data = {'Algorithm': [], 'TumorType': [], 'TestAUC': []}

for cancer, auc_dict in zip(cancer_types, auc_values):
    for algo, auc in auc_dict.items():
        data['Algorithm'].append(algo)
        data['TumorType'].append(cancer)
        data['TestAUC'].append(auc)

df = pd.DataFrame(data)

# Plot the boxplot using Seaborn
plt.figure(figsize=(12, 6))
sns.boxplot(x='Algorithm', y='TestAUC', data=df, color="skyblue")

# Add points and labels for each cancer type
sns.swarmplot(x='Algorithm', y='TestAUC', data=df, hue='TumorType', size=8, edgecolor="gray", linewidth=1)
plt.legend(title='Tumor Type', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.xlabel('Algorithm')
plt.ylabel('Test AUC')
plt.title('Test AUC by Algorithm and Tumor Type')

# Add the name of the data used
data_name = "YOUNG (MC38), TOP 100 & BY CANCER TYPE"
plt.text(0, -0.15, f"Data: {data_name}", transform=plt.gca().transAxes, fontsize=12)

plt.show()
