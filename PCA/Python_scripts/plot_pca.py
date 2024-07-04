import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import sys

pca_data = pd.read_csv(sys.argv[1], sep=',')
pca_country = pd.read_csv(sys.argv[2])
pca_country = pd.read_csv(sys.argv[2], sep='\t')
"""
pca_country['country'] = pca_country['country'].map(lambda x: 'Central Africa' if x in ['Gabon', 'Cameroon', 'Democratic Republic of the Congo',] else 
                                                    ('West Africa' if x in ['Ghana', 'Nigeria', 'Mauritania', 'Guinea', "Côte d'Ivoire", "Mali", 'Burkina Faso', 'Senegal', 'Gambia', 'Benin'] else 
                                                     ('East Africa' if x in ['Ethiopia', 'Tanzania', 'Kenya', 'Mozambique', 'Madagascar', 'Uganda', 'Sudan', 'Malawi' ] else ('Sao Tome' if x in ['Sao_Tome'] else print()))))
"""
# Create a scatter plot
country = pca_country[['sample', 'country']]
df_pca = pd.merge(pca_data, country, on="sample")

# Get unique clusters
unique_clusters = df_pca['country'].unique()

# Get the tab10 palette colors
tab10_colors = sns.color_palette('tab10')

# Define colors for each group using tab10 colors
color_dict = {'Kenya': tab10_colors[0], 'Gabon': "#FFA500", 'Ghana': '#B2D63F', 'Sao_Tome': "#FF0000", 'Tanzania': tab10_colors[4], 'Cameroon': "green", 'Nigeria': "#00D000"}


# Set a color palette with a unique color for each cluster
cluster_palette = sns.color_palette('tab10', n_colors=len(unique_clusters))

# Define your own color palette
custom_colors = ['yellow', 'purple', 'teal', 'red']

plt.figure(figsize=(9, 5))
sns.scatterplot(x='PC1', y='PC2', data=df_pca, hue='country', palette=color_dict, s=30)
plt.xlabel('PC1 (6.08%)', fontsize=16, labelpad=10, fontweight='bold')
plt.ylabel('PC2 (3.74%)', fontsize=16, labelpad=10, fontweight='bold')
legend = plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False, fontsize=14)  # Adjust legend position
plt.gca().tick_params(axis='both', which='major', width=2)
plt.xticks( fontsize=14, fontweight='bold')
plt.yticks( fontsize=14, fontweight='bold')
plt.gca().spines['bottom'].set_color('black')  # Adjust bottom axis line color
plt.gca().spines['bottom'].set_linewidth(2)    # Adjust bottom axis line width
plt.gca().spines['left'].set_color('black')    # Adjust left axis line color
plt.gca().spines['left'].set_linewidth(2)
# Remove the frame around the plot
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.tight_layout()

# Modify legend labels
for text in legend.get_texts():
    if text.get_text() == 'Sao_Tome':
        text.set_text('STP')
plt.show()
