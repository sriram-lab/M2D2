# GRAPHS OF DATA 
data_s2 = pd.read_excel("data_s2.xlsx")
ma = pd.read_excel("data_Ma.xlsx")
rhoads = pd.read_excel("data_Rhoads.xlsx")
yilancioglu = pd.read_excel("data_Yilancioglu.xlsx")

# Histogram of all interaction scores
plt.hist(data_s2['score'])
plt.title('Histogarm of All Interaction Scores (data_s2.xlsx, all literatures)')
plt.xlabel('Interaction Score')
plt.ylabel('Count')

# Histogram of range of scores for each literature/data set
data_s2['Dataset'] = 'All'
ma['Dataset'] = 'Ma'
rhoads['Dataset'] = 'Rhoads'
yilancioglu['Dataset'] = 'Yilancioglu'

combined = pd.concat([
    data_s2[['score', 'Dataset']],
    ma[['score', 'Dataset']],
    rhoads[['score', 'Dataset']],
    yilancioglu[['score', 'Dataset']]
])

plt.figure(figsize=(8, 5))
plt.title('Interaction Scores Across Datasets')
plt.hist(combined[combined['Dataset'] == 'All']['score'], alpha=0.6, label='All datasets (data_s2.xlsx) (332)')
plt.hist(combined[combined['Dataset'] == 'Ma']['score'], alpha=0.6, label='Ma 2019 (241)')
plt.hist(combined[combined['Dataset'] == 'Rhoads']['score'], alpha=0.6, label='Rhoads (39)')
plt.hist(combined[combined['Dataset'] == 'Yilancioglu']['score'], alpha=0.6, label='Yilancioglu 2019 (120)')
plt.xlabel('Interaction Score')
plt.ylabel('Count')
plt.legend()
plt.show()

plt.figure(figsize=(8, 5))
plt.title('Interaction Scores Across Datasets')
plt.hist(combined[combined['Dataset'] == 'All']['score'], alpha=0.6, label='All datasets (data_s2.xlsx) (332)', density=True)
plt.hist(combined[combined['Dataset'] == 'Ma']['score'], alpha=0.6, label='Ma 2019 (241)', density=True)
plt.hist(combined[combined['Dataset'] == 'Rhoads']['score'], alpha=0.6, label='Rhoads (39)', density=True)
plt.hist(combined[combined['Dataset'] == 'Yilancioglu']['score'], alpha=0.6, label='Yilancioglu 2019 (120)', density=True)
plt.xlabel('Interaction Score')
plt.ylabel('Probability Density') # Because the datasets are different sizes, like 39 in Rhoads vs 241 in Ma
plt.legend()
plt.show()

# Boxplot of interaction score ranges for the repeats
plt.figure(figsize=(8, 5))
data_s2['combinations'] = data_s2['Drug_1'] + ' + ' + data_s2['Drug_2']
data_s2.boxplot(column='score', by='combinations', rot=45)
plt.title('Boxplot of Interaction Score Ranges from each Drug Pair (data_s2.xlsx, all literatures)')
plt.xlabel('Drug Combination')
plt.ylabel('Interaction Score')
plt.show()

# Heatmaps for synergies
# cool colors: rocket_r, coolwarm
plt.figure(figsize=(8, 5))
pivot_s2 = data_s2.pivot_table(values='score', index='Drug_1', columns='Drug_2')
sns.heatmap(pivot_s2, cmap='coolwarm')
plt.title('data_s2.xlsx')
plt.show()

plt.figure(figsize=(8, 5))
pivot_rhoads = rhoads.pivot_table(values='score', index='Drug_1', columns='Drug_2')
sns.heatmap(pivot_rhoads, cmap='coolwarm')
plt.title('Rhoads (39)')
plt.show()

plt.figure(figsize=(8, 5))
pivot_ma = ma.pivot_table(values='score', index='Drug_1', columns='Drug_2')
sns.heatmap(pivot_ma, cmap='coolwarm')
plt.title('Ma 2019 (241)')
plt.show()

plt.figure(figsize=(8, 5))
pivot_yilancioglu = yilancioglu.pivot_table(values='score', index='Drug_1', columns='Drug_2')
sns.heatmap(pivot_yilancioglu, cmap='coolwarm')
plt.title('Yilancioglu 2019 (120)')
plt.show()
