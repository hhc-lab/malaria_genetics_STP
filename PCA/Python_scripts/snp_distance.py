import sys
import numpy as np
import pandas as pd
from Bio import AlignIO
from Bio.Seq import Seq
import subprocess
from Bio.Align import MultipleSeqAlignment 


def calculate_average_snp_distance(seqs):
    num_seqs = len(seqs)
    total_distance = 0
    pair_count = 0
    
    # Iterate over every pair of sequences
    for i in range(num_seqs):
        for j in range(i+1, num_seqs):  # Avoid redundant comparisons
            distance = snp_distance(seqs[i], seqs[j])
            total_distance += distance
            pair_count += 1
    
    # Calculate the average SNP distance
    if pair_count > 0:
        avg_distance = total_distance / pair_count
    else:
        avg_distance = 0
    
    return avg_distance


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
    filtered_sequences = []
    for record in alignment:
        filtered_seq = ''
        for i in range(len(homogeneity)):
            if  gap_percentage[i] <= 0.05:
#            if homogeneity[i] <= 0.8 and gap_percentage[i] <= 0.5:
                filtered_seq += record.seq[i]
        filtered_record = record[:]
        filtered_record.seq = Seq(filtered_seq)
        filtered_sequences.append(filtered_record)
    filtered_alignment = MultipleSeqAlignment(filtered_sequences)
    return filtered_alignment


# Define the snp_distance function
def snp_distance(seq1, seq2, gap_char='-'):
    """
    Calculate the SNP distance between two sequences, skipping positions with gaps.

    Parameters:
    seq1 (str): The first sequence.
    seq2 (str): The second sequence.
    gap_char (str): The character representing a gap. Default is '-'.

    Returns:
    int: The SNP distance between the two sequences, excluding gaps.
    
    Raises:
    ValueError: If the sequences are not of the same length.
    """
    
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of the same length")

    # Initialize the SNP distance counter
    distance = 0
    position = 0
    # Iterate through the sequences and count the differences, skipping gaps
    for base1, base2 in zip(seq1, seq2):
        if base1 == gap_char or base2 == gap_char:
            continue
        position += 1
        if base1 != base2:
            distance += 1
    distance_percentage = distance/position

    print(seq1, seq2, "SNP distance: ",distance)
    return distance_percentage


# Function to read the FASTA alignment file and convert it to a pandas DataFrame
def alignment_to_dataframe(file_path):
    alignment = AlignIO.read(file_path, "fasta")
    data = []
    for record in alignment:
        sample_name = record.id
        sequence = str(record.seq)
        data.append([sample_name, sequence])
    df = pd.DataFrame(data, columns=["sample", "sequence"])
    return df

def calculate_country_snp_distances(sample_id, country, data):
    country_data=data[data['country'] == country]
    sample_seq = data.loc[data['sample']== sample_id, 'sequence'].iloc[0]
    sum_distance = 0
    count_pairs = 0
    output = pd.DataFrame(columns=["sequence1", "sequence2", "snp_distance"])    
    for _, row in country_data.iterrows():
        if row['sample'] != sample_id:
            seq2= row['sequence']
            distance = snp_distance(sample_seq, seq2)
            sum_distance += distance
            count_pairs += 1
            new_record = pd.DataFrame([{"sequence1":sample_id, "sequence2":row['sample'], "snp_distance":distance}])
            output = pd.concat([output, new_record], ignore_index=True)
    if count_pairs > 0:
        average_snp_distance= sum_distance/count_pairs
        file_name = f"Results/snp_distance/{sample_id}_{country}.csv"
        output.to_csv(file_name, index=False)
    else:
        print("Error no sequences found for ", country)
    return average_snp_distance 

# Read the sequences and store them into a pandas dataframe
data = pd.read_csv(sys.argv[1], sep='\t')
print("data read from file",data)

# Write sequences to a temporary FASTA file
with open("sequences.fasta", "w") as fasta_file:
    for i, row in data.iterrows():
        fasta_file.write(f">{row['sample']}\n{row['sequence']}\n")

# Using multseq malign
print('Performing multyple sequence alignment')
multseq_cmd = f"multseq malign --in sequences.fasta --threads 12 --package muscle --out alignment.fasta --verbose"
subprocess.run(multseq_cmd, shell=True, check=True)

# Read the data
output_file = "filtered_alignment.fasta"

# Read alignment
alignment = AlignIO.read("alignment.fasta", "fasta")

# Filter alignment
filtered_alignment = filter_alignment(alignment)

# Write filtered alignment to output file
AlignIO.write(filtered_alignment, output_file, "fasta")

print("Filtered alignment saved to", output_file)

filtered_alignment = alignment_to_dataframe("filtered_alignment.fasta")
data.drop(columns=["sequence"], inplace=True)
alignment_data= pd.merge(filtered_alignment, data, on="sample", how="left",)
sao_tome_data = alignment_data[data['country'] == 'Sao_Tome']
print('sao_tome_data:',sao_tome_data)
# Initialize a dictionary to store average SNP distances for each country
country_distances = pd.DataFrame(columns=['sample', 'country', 'snp_distance']) 
countries = alignment_data['country'].unique()
# Calculate average SNP distance for each sample in Sao_Tome against all other samples in the same country
for index, row in sao_tome_data.iterrows():
    sample_id = row['sample']
    for c in countries:
        average_distance = calculate_country_snp_distances(sample_id, c, alignment_data)
        new_record = pd.DataFrame([{'sample': sample_id, 'country': c, 'snp_distance': average_distance}])
        country_distances = pd.concat([country_distances, new_record], ignore_index=True)
stp_seqs = sao_tome_data['sequence']
seqs = stp_seqs.tolist()
stp_average_distance = calculate_average_snp_distance(seqs)
print("STP average SNP distance:",stp_average_distance)

# Save the table to a CSV file
country_distances.to_csv('country_average_snp_distances.csv', index=False)

