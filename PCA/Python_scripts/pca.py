import sys
import numpy as np
from sklearn.decomposition import PCA
import pandas as pd
from Bio import AlignIO
from Bio.Seq import Seq
import sys
import subprocess
from Bio.Align import MultipleSeqAlignment 
from Bio.SeqRecord import SeqRecord



def one_hot_encode(sequence):
    nucleotides = {'A': 0, 'T': 1, 'C': 2, 'G': 3}
    encoding = np.zeros((len(sequence), len(nucleotides) + 1))  # Additional index for other symbols
    other_symbols_count = 0  # Counter for other symbols
    other_symbols_percentage = []

    for i, nucleotide in enumerate(sequence):
        upper_nucleotide = nucleotide.upper()
        if upper_nucleotide in nucleotides:
            encoding[i, nucleotides[upper_nucleotide]] = 1
        else:
            encoding[i, -1] = 1  # Encoding other symbols as 4
            other_symbols_count += 1  # Increment the counter for other symbols
    
    # Calculate the percentage of other symbols
        percentage_other_symbols = (other_symbols_count / len(sequence)) * 100
    
    return encoding.flatten(), percentage_other_symbols

def read_sequence_table(file_path):
    df = pd.read_table(file_path, sep='\t')
    return df

def calculate_homogeneity(alignment):
    homogeneity = []
    gap_percentage = []
    alignment_length = alignment.get_alignment_length()
    for i in range(alignment_length):
        column = alignment[:, i]
        gap_count = column.count('-')
        gap_percentage.append(gap_count / len(column))
        counts = {x: column.count(x) for x in set(column)}
        max_count = max(counts.values())
        homogeneity.append(max_count / len(column))
    print("homogeneity:", homogeneity)
    print("gap_percentage:", gap_percentage)
    return homogeneity, gap_percentage

def filter_alignment(alignment):
    homogeneity, gap_percentage = calculate_homogeneity(alignment)
    filtered_alignment = MultipleSeqAlignment([])
    for i in range(len(alignment)):
        filtered_seq = ''
        for j in range(len(homogeneity)):
            if homogeneity[j] <= 0.8 and gap_percentage[j] <=0.001:
                filtered_seq += alignment[i][j]
        print(filtered_seq)
        filtered_seq_record = SeqRecord(Seq(filtered_seq), id=alignment[i].id, description=alignment[i].description)
        filtered_alignment.append(filtered_seq_record)

    return filtered_alignment

def random_sample(df, n):
    """
    Takes a random sample of each country from the DataFrame.
    
    Returns:
    DataFrame: A DataFrame containing the random samples for each country.
    """
    # Function to sample each group
    def sample_group(group):
        return group.sample(min(len(group), n))
    
    # Group by 'country' and apply the sampling function to each group
    sampled_df = df.groupby('country').apply(sample_group).reset_index(drop=True)
    
    return sampled_df


# Read the sequences and store them into a pandas dataframe
raw_data= read_sequence_table(sys.argv[1])
print("data read from file", raw_data)

sample_no = int(sys.argv[2])
data = random_sample(raw_data, sample_no)


# Write sequences to a temporary FASTA file
with open("sequences.fasta", "w") as fasta_file:
    for i, row in data.iterrows():
        fasta_file.write(f">{row['sample']}\n{row['sequence']}\n")


# Using multseq malign
print('Performing multyple sequence alignment')
multseq_cmd = f"multseq malign --in sequences.fasta --threads 12 --package muscle --out alignment.fasta --verbose"
subprocess.run(multseq_cmd, shell=True, check=True)


output_file = "filtered_alignment.fasta"

# Read alignment
alignment = AlignIO.read("alignment.fasta", "fasta")

# Filter alignment
filtered_alignment = filter_alignment(alignment)

# Write filtered alignment to output file
AlignIO.write(filtered_alignment, output_file, "fasta")

print("Filtered alignment saved to", output_file)

filtered_alignment = AlignIO.read("filtered_alignment.fasta", "fasta")

# Convert alignment to a list of sequences
sequences = [str(record.seq) for record in filtered_alignment]
print(sequences)

# One-Hot Encoding the sequences
print('One-hot encodding the sequences')
#X = np.array([one_hot_encode(seq) for seq in sequences])

# Process the sequences and capture the encoded sequences and percentages
encoded_sequences = []
other_symbols_percentage = []

for seq in sequences:
    encoded_seq, percentage = one_hot_encode(seq)
    encoded_sequences.append(encoded_seq)
    other_symbols_percentage.append(percentage)

# Convert the encoded sequences to a NumPy array
X = np.array(encoded_sequences)


# Calculate the average percentage of other symbols over all sequences
average_percentage = sum(other_symbols_percentage) / len(other_symbols_percentage)
print("Average percentage of other symbols in the sequences:",average_percentage)


other_symbols_percentage = pd.DataFrame(other_symbols_percentage)
other_symbols_percentage.to_csv("other_symbols_percentage.csv", index=False)



# PCA
print('Performing PCA')
pca = PCA(n_components=10)
X_pca = pca.fit_transform(X)

# Create a DataFrame with the PCA-transformed data
df_pca = pd.DataFrame(X_pca, columns=[f'PC{i+1}' for i in range(X_pca.shape[1])])
df_pca['sample'] = data['sample']

# Add explained variance ratio to DataFrame
explained_variance_ratio = pca.explained_variance_ratio_

# Save to CSV
df_variance=pd.DataFrame(explained_variance_ratio)
df_variance.to_csv('variance.csv', index=False)

# Save to CSV
df_pca.to_csv('pca_transformed_data.csv', index=False)

# Save samples to file
data.to_csv('pca_input_data.csv', index=False)
