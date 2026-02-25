import tkinter as tk
seq='''GTGACCCTGGCCAGGACTGACCTGGAGATGCAGATCGAAGGCCTGAAGGAGGAGCTGGCC
TACCTGAGGAAGAACCACGAGGAGGAGATGCTTGCTCTGAGAGGTCAGACCGGCGGAGAT
GTGAACGTGGAGATGGATGCTGCACCTGGCGTGGACCTGAGCCGCATCCTGAATGAGATG
CGTGACCAGTACGAGCAGATGGCAGAGAAAAACCGCAGAGACGCTGAGACCTGGTTCCTG
AGCAAGACCGAGGAGCTGAACAAAGAAGTGGCCTCCAACAGCGAACTGGTACAGAGCAGC
CGCAGTGAGGTGACGGAGCTCCGGAGGGTGCTCCAGGGCCTGGAGATTGAGCTGCAGTCC
CAGCTCAGCACGAAAGCATCCCTGGAGAACAGCCTGGAGGAGACCAAAGGCCGCTACTGC
ATGCAGCTGTCCCAGATCCAGGGACTGATTGGCAGTGTGGAGGAGCAGCTGGCCCAGCTA
CGCTGTGAGATGGAGCAGCAGAGCCAGGAGTACCAGATCTTGCTGGATGTGAAGACGCGG
CTGGAGCATGAGATTGCCACCTACCGCCGCCTGCTGGANGGCGAGGATGCCCACCTTTCC
TCCCAGCAAGCATCTGGCCAATCCTATTCTTCCCGCGAGGTCTTCACCTCCTCCTCGTCC
TCTTCGAGCCGTCAGACCCGACCCATCCTCAAGGAGCAGAGCTCATCCAGCTTCAGCCAG
GGCCAGAGTTCCTAG'''
#FOR SEARCHING U AND A IN THE SEQUENCE IF ITS DNA OR RNA
#SEQUENCE IDENTIFICATION
def search(seq):
    if "T" in seq:
        print("ITS DNA SEQUENCE")
    elif "U" in seq:
        print("ITS RNA SEQUENCE")
    else:
        print("INVALID SEQUENCE")
search(seq)
 # THIS PART IS FOR TRANSCRIBING THE SEQUENCE BY USING FUNCTION  
#TRANSCRIPTION 
def replace_bases(seq):
    dict = {'T': 'U', 'A': 'T', 'C': 'G', 'G': 'C'}

    rna_sequence = ''
    for base in seq:
        if base in dict:
           rna_sequence += dict[base]
        else:
            rna_sequence += base  # If base is not A, T, C, or G, keep it unchanged

    return rna_sequence 
RNA_sequence = replace_bases(seq)
print("THE GIVEN SEQUENCE IS NOW RNA SEQUENCE:", RNA_sequence)
#CODONS OF THE INPUT SEQUENCE( WHICH WILL HELP IN MAKING THE ORFS(ON BASES OF START AND STOP CODONS))
def extract_codons(seq):
    codons = [seq[i:i+3]
            for i in range(0, len(seq), 3)]
    return codons
codons = extract_codons(seq)
print("THE DNA SEQUENCE IS CONVERTED INTO CODONS", codons)
#ORFS OF THE SEQUENCE
def find_orfs(seq):
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    orfs = []

    # Finding start codons
    for i in range(len(seq) - 2):
        if seq[i:i + 3] == start_codon:
            # Finding stop codons after start codon
            for j in range(i + 3, len(seq), 3):
                codon = seq[j:j + 3]
                if codon in stop_codons:
                    orfs.append(seq[i:j + 3])
                    break

    return orfs

# Find ORFs
resulting_orfs = find_orfs(seq)

# Print the found ORFs
for index, orf in enumerate(resulting_orfs, 1):
    print(f"ORF {index}: {orf}")
#TRANSLATION (the process in which the transcribe sequence making the protein sequence)
def translate_sequence(RNA_sequence):
    codon_table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '', 'TAG': '',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    }

    protein_sequence = ''
    for i in range(0, len(RNA_sequence), 3):
        codon = RNA_sequence[i:i + 3]
        amino_acid = codon_table.get(codon, '')
        if amino_acid == '*':
            break  # Stop codon encountered, terminate translation
        protein_sequence += amino_acid

    return protein_sequence

# Translate the DNA sequence into amino acids
resulting_protein_sequence = translate_sequence(RNA_sequence)

print(f"Translated protein sequence: {resulting_protein_sequence}")
  
def process_sequence():
    seq = sequence_entry.get()

    result_text.delete(1.0, tk.END)  # Clear previous content
    result_text.insert(tk.END, f"Processed sequence: {seq}\n\n")
    #result_text.insert(tk.END, f"{search(seq)}\n\n")
    result_text.insert(tk.END, f"THE SEQUENCE IS {'RNA' if 'U' in seq else 'DNA' if 'T' in seq else 'INVALID'}\n\n")
   
    RNA_sequence = replace_bases(seq)
    result_text.insert(tk.END, f"The given sequence is now RNA sequence: {RNA_sequence}\n\n")
    
    codons = extract_codons(seq)
    result_text.insert(tk.END, f"The DNA sequence is converted into codons: {codons}\n\n")
    
    orfs = find_orfs(seq)
    result_text.insert(tk.END, f"\nOpen Reading Frames (ORFs):\n")
    for index, orf in enumerate(orfs, 1):
        result_text.insert(tk.END, f"ORF {index}: {orf}\n\n")
    
    translated_sequence = translate_sequence(RNA_sequence)
    result_text.insert(tk.END, f"\nTranslated Protein Sequence:\n{translated_sequence}")

root = tk.Tk()
root.title("Sequence Analyzer")

input_frame = tk.Frame(root)
input_frame.pack(padx=10, pady=10)

sequence_label = tk.Label(input_frame, text="Enter Sequence:")
sequence_label.pack()

sequence_entry = tk.Entry(input_frame, width=60)
sequence_entry.pack()

process_button = tk.Button(root, text="Process Sequence", command=process_sequence)
process_button.pack(pady=10)

result_text = tk.Text(root, height=20, width=80)
result_text.pack(padx=10, pady=10)

root.mainloop()