import sys
import pandas as pd

# Read the input CSV file
country_distances = pd.read_csv(sys.argv[1])

# Get unique sample IDs
unique_sample_ids = country_distances['sample'].unique()

# Initialize an empty list to store the results
result = []

# Iterate over each unique sample ID
for sample_id in unique_sample_ids:
    # Find the SNP distance for Sao Tome
    san_tome_snp_distance = country_distances.loc[(country_distances['sample'] == sample_id) & (country_distances['country'] == 'Sao_Tome'), 'snp_distance'].values[0]
    
    # Find the SNP distances for all other countries
    other_countries_data = country_distances.loc[(country_distances['sample'] == sample_id) & (country_distances['country'] != 'Sao_Tome')]
    
    # Store the result for each sample ID, country, SNP distance, and other columns
    for index, row in other_countries_data.iterrows():
        if row['snp_distance'] < san_tome_snp_distance:
            result.append({'sample': sample_id,
                           'country': row['country'],
                           'snp_distance': row['snp_distance'],
                           'Sao_Tome_snp_distance': san_tome_snp_distance,
                           'snp_distance_smaller_than_sao_tome': True})
        else:
            result.append({'sample': sample_id,
                           'country': row['country'],
                           'snp_distance': row['snp_distance'],
                           'Sao_Tome_snp_distance': san_tome_snp_distance,
                           'snp_distance_smaller_than_sao_tome': False})

# Convert the result list to a DataFrame
result_df = pd.DataFrame(result)

# Print the result DataFrame
print(result_df)

# Save the result DataFrame to a CSV file
result_df.to_csv('imported.csv', index=False)



"""import sys
import pandas as pd

country_distances = pd.read_csv(sys.argv[1])



unique_sample_ids = country_distances['sample'].unique()

# Initialize an empty list to store the resultsi
result = []

# Iterate over each unique sample ID
for sample_id in unique_sample_ids:
    # Find the SNP distance for San Tome
    san_tome_snp_distance = country_distances.loc[(country_distances['sample'] == sample_id) & (country_distances['country'] == 'Sao_Tome'), 'snp_distance'].values[0]
    
    # Find the SNP distances for all other countries
    other_countries_snp_distances = country_distances.loc[(country_distances['sample'] == sample_id) & (country_distances['country'] != 'Sao_Tome'), 'snp_distance']
    
    # Check if any other country has a smaller SNP distance than San Tome
    if (other_countries_snp_distances < san_tome_snp_distance).any():
        result.append({'sample': sample_id, 'snp_distance_smaller_than_san_tome': True})
    else:
        result.append({'sample': sample_id, 'snp_distance_smaller_than_san_tome': False})

# Convert the result list to a DataFrame
result_df = pd.DataFrame(result)

# Print the result DataFrame
print(result_df)


result_df.to_csv('imported.csv', index=False)
"""
