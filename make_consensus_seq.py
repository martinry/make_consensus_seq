seqs = []

with open('msa.txt') as sequences:
    for s in sequences:
        s = s.strip('\n')
        if(not s.startswith('>')): # skip lines with '>' in case of fasta
            seqs.append(s)

letters = ['A', 'C', 'G', 'T']

iupac_codes = {
     'A': [1, 0, 0, 0],
     'C': [0, 1, 0, 0],
     'G': [0, 0, 1, 0],
     'T': [0, 0, 0, 1],
     'W': [1, 0, 0, 1],
     'S': [0, 1, 1, 0],
     'Y': [0, 1, 0, 1],
     'R': [1, 0, 1, 0],
     'M': [1, 1, 0, 0],
     'K': [0, 0, 1, 1],
     'H': [1, 1, 0, 1],
     'B': [0, 1, 1, 1],
     'V': [1, 1, 1, 0],
     'D': [1, 0, 1, 1],
     'N': [1, 1, 1, 1]
     }

seq_length = len(seqs[0])

transposed_seqs = []


for i in range(0,seq_length): # Count from zero to the number of nucleotides in the seq
    nucleotides = [n[i] for n in seqs]  # For every sequence, take the nucleotide for the
                                        # current iteration
    transposed_seqs.append(nucleotides) # We now have a transposed list of the sequences
                                        # (switch rows and columns)
    
summary_table = []  # This will contain a binary representation of what each column contains
                    # in accordance to the iupac_codes dictionary

for subseq in transposed_seqs:  # Iterate through each seq (ie column in the original dataset)
    temp = [0,0,0,0]               # Empty binary array for each iteration
    seqstr = ''.join(subseq)    # Join inner elements of nested list
    for i in letters:           # For every possible nucleotide
        ind = letters.index(i)  # find the index of that nucleotide
        if(i in seqstr):        # if the current nucl of 'letters' exists in our sequence
            temp[ind] = 1       # then set that index's value to 1
    summary_table.append(temp)     # and append it to our binary table

consensus = [] # The full consensus string

for i in summary_table:
    for k,v in iupac_codes.items():
        if i == v:
            consensus.extend(k)

print(''.join(consensus))
