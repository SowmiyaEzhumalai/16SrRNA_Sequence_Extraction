import os
from Bio import SeqIO

def extract_16S_rRNA_from_folder(input_folder, output_folder):
    # Ensure output folder exists
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    # Iterate over all GenBank files in the input folder
    for file_name in os.listdir(input_folder):
        if file_name.endswith(".gb") or file_name.endswith(".gbk"):  # Check for GenBank file extensions
            input_file = os.path.join(input_folder, file_name)
            with open(input_file, "r") as gb_file:
                for record in SeqIO.parse(gb_file, "genbank"):
                    # Try to extract strain name from the GenBank metadata
                    strain_name = record.annotations.get("organism", "unknown_strain").replace(" ", "_")
                    
                    # Search for 16S rRNA feature
                    for feature in record.features:
                        if feature.type == "rRNA" and "16S" in feature.qualifiers.get("product", [""])[0]:
                            # Extract the 16S rRNA sequence
                            rRNA_seq = feature.location.extract(record.seq)
                            #Saving input file name without the extension in a variable named 'file_base_name'
                            file_base_name = os.path.splitext(file_name)[0]
                            # Create output file path using strain name
                            output_file = os.path.join(output_folder, f"{file_base_name}.fasta")
                            
                            # Write the sequence to the output file in FASTA format
                            with open(output_file, "w") as output_fasta:
                                output_fasta.write(f">{file_base_name}\n")
                                output_fasta.write(str(rRNA_seq) + "\n")
                            print(f"16S rRNA sequence extracted and saved to {output_file}")
                            break  # Stop after finding the first 16S rRNA feature
                    else:
                        print(f"No 16S rRNA sequence found in {file_name}")

# Example usage
input_genbank_folder = "Path/to/your/input/directory/"
output_fasta_folder = "Path/to/your/output/directory/"
extract_16S_rRNA_from_folder(input_genbank_folder, output_fasta_folder)
