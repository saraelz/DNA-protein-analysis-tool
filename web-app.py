
######################
# Import libraries
######################

import pandas as pd
import streamlit as st
import altair as alt
from PIL import Image

######################
# Page Title
######################

image = Image.open('central-dogma.png')

st.image(image, use_column_width=True)

st.write("""
# DNA to Protein

This app analyzes the nucleotide composition of query DNA.
         
Using the central dogma of biology, we will convert DNA to mRNA to Protein

***
""")


######################
# Input Text Box
######################

#st.sidebar.header('Enter DNA sequence')
st.header('Enter DNA sequence')

sequence_input = "TACGAACACGTGGAGGCAAACAGGAAGGTGAAGAAGAACTTATCCTATCAGGACGGAAGGTCCTGTGCTCGGG\nATCTTCCAGACGTCGCGACTCTAAATTGCCCCCTCTGAGGTCAAGGAACACAAGATGGTTTTGGAAATGC\nTGAACCCGATACATTATAACATCACCAGCATCGTGCCTGAAGCCATGCCTGCTGCCACCATGCCAGTCCT"

#sequence = st.sidebar.text_area("Sequence input", sequence_input, height=250)
user_input = st.text_area("Please input DNA sequence from 3' to 5'", sequence_input, height=250)


def check_DNA_validity(input):
  clean_input = []
  for char in input:
    if char.isalpha():
      char = char.upper()
      if char in ['A', 'T', 'G', 'C']:
        clean_input.append(char)
      else:
        # Check DNA
        st.warning('Warning: Input DNA should only include the letters A, T, G, and C.')
  return ''.join(clean_input)    

sequence = check_DNA_validity(user_input)

st.write("""
***
""")

## Prints the input DNA sequence
st.header('INPUT (DNA Query)')
st.write('3\' ' + sequence + ' 5\' ')

## DNA nucleotide count
st.header('OUTPUT (DNA Nucleotide Count)')

### 1. Print dictionary
st.subheader('1. Print dictionary')
def DNA_nucleotide_count(seq):
  d = dict([
            ('A',seq.count('A')),
            ('T',seq.count('T')),
            ('G',seq.count('G')),
            ('C',seq.count('C'))
            ])
  return d

X = DNA_nucleotide_count(sequence)

X

### 2. Print text
st.subheader('2. Print text')
st.write('There are  ' + str(X['A']) + ' adenine (A)')
st.write('There are  ' + str(X['T']) + ' thymine (T)')
st.write('There are  ' + str(X['G']) + ' guanine (G)')
st.write('There are  ' + str(X['C']) + ' cytosine (C)')
st.write(f'There are {len(sequence)} total nucleotides.')

### 3. Display DataFrame
st.subheader('3. Display DataFrame')
df = pd.DataFrame.from_dict(X, orient='index')
df = df.rename({0: 'count'}, axis='columns')
df.reset_index(inplace=True)
df = df.rename(columns = {'index':'nucleotide'})
st.write(df)

### 4. Display Bar Chart using Altair
st.subheader('4. Display Bar chart')
p = alt.Chart(df).mark_bar().encode(
    x='nucleotide',
    y='count'
)
p = p.properties(
    width=alt.Step(80)  # controls width of bar.
)
st.write(p)


### 5. Transcription and Translation
st.subheader('5. DNA to mRNA to Protein')

st.write('The input DNA query is used as a template for transcription. Then, mRNA is used to encode proteins. Keep in mind that adenine binds to uracil in RNA.')
'Input DNA' 
st.code('3\'' + sequence + '5\' ')

# Convert DNA to mRNA
transcription_mappings = {
  'A': 'U',
  'T': 'A',
  'G': 'C',
  'C': 'G'
}
mRNA_ribonucleotides = []
for nucleotide in sequence:
  if nucleotide in transcription_mappings:
    mRNA_ribonucleotides.append(transcription_mappings[nucleotide])
mRNA_sequence = ''.join(mRNA_ribonucleotides)
'mRNA'
st.code('5\'' + mRNA_sequence + '3\' ')

# Convert DNA to tRNA
tRNA_sequence = sequence.replace('T', 'U')
'tRNA'
st.code('3\'' + tRNA_sequence + '5\' ')

# Convert mRNA to Protein (translaltion)
RNA_codon_table = {
# U
'UUU': 'Phe', 'UCU': 'Ser', 'UAU': 'Tyr', 'UGU': 'Cys', # UxU
'UUC': 'Phe', 'UCC': 'Ser', 'UAC': 'Tyr', 'UGC': 'Cys', # UxC
'UUA': 'Leu', 'UCA': 'Ser', 'UAA': 'STOP', 'UGA': 'STOP', # UxA
'UUG': 'Leu', 'UCG': 'Ser', 'UAG': 'STOP', 'UGG': 'Trp', # UxG

# C
'CUU': 'Leu', 'CCU': 'Pro', 'CAU': 'His', 'CGU': 'Arg', # CxU
'CUC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg', # CxC
'CUA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg', # CxA
'CUG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg', # CxG

# A
'AUU': 'Ile', 'ACU': 'Thr', 'AAU': 'Asn', 'AGU': 'Ser', # AxU
'AUC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser', # AxC
'AUA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg', # AxA
'AUG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg', # AxG

# G
'GUU': 'Val', 'GCU': 'Ala', 'GAU': 'Asp', 'GGU': 'Gly', # GxU
'GUC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly', # GxC
'GUA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly', # GxA
'GUG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'  # GxG
}


singleletter = {'Cys': 'C', 'Asp': 'D', 'Ser': 'S', 'Gln': 'Q', 'Lys': 'K',
'Trp': 'W', 'Asn': 'N', 'Pro': 'P', 'Thr': 'T', 'Phe': 'F', 'Ala': 'A',
'Gly': 'G', 'Ile': 'I', 'Leu': 'L', 'His': 'H', 'Arg': 'R', 'Met': 'M',
'Val': 'V', 'Glu': 'E', 'Tyr': 'Y', 'STOP': '*'}

def find_str_in_Str(str, start_codon):
  """ 
    Find the index of the first occurrence in a String,

    str: string of ribonucleotides (mRNA) containing only A,U,C,G
    start_codon: list of strings
    
    Usage: Translation will begin at the index of the first start codon in the mRNA
  """
  slow = 0
  while slow < len(str):
    fast = slow
    for i in range(3):
      fast = slow + i
      if start_codon[i] != str[fast]:
        break
      if i==2 and start_codon[i] == str[fast]:
        return slow
    slow = slow + 1
  return -1
    

def translation(mRNA, start_index):
  # Assume mRNA is a string with only characters A U C G

  codons = [mRNA[i:i+3] for i in range(start_index,len(sequence), 3)]
  'Codons'
  st.code(codons)

  amino_acids = []
  for codon in codons:
    if codon in RNA_codon_table:
      amino_acids.append(RNA_codon_table[codon])
    else:
      amino_acids.append(codon)
  return amino_acids


start_index = find_str_in_Str(mRNA_sequence, 'AUG')
amino_acids = translation(mRNA_sequence,start_index)

'Amino Acids'
st.code(amino_acids)

# Get single letter amino acids
amino_acids_sl = [singleletter[aa] for aa in amino_acids if singleletter.get(aa)]
amino_acids_sl_str = ''.join(amino_acids_sl)
'Single Letter Amino Acids'
st.code(amino_acids_sl_str)


'Stop codons do not code for any amino acids, and they are represented by the \'*\'. For the purposes of this learning module, the stop codons have been ignored.'

'Credits'
'This program was written by Sara Elzeiny. This program is intended for learning purposes only. It is not intended for professional use.'