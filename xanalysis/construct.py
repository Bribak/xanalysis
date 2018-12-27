from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import pandas as pd
import numpy as np

# function to transform cDNA library into protein library
def translate_aa(filename, outputfile):
    with open(filename) as inputfile:
        for line in inputfile:
                line = line.replace(' ','')
                line = line.replace('\n','')
                line = Seq(line, IUPAC.unambiguous_dna)
                line = str(line.translate())
                line = line.replace('*','')
                with open(outputfile,'a') as the_file:
                    the_file.write(line+'\n')

# function to delete fasta ids, short ORFs and splits dna sequences into triplets
# one gene per line
def prep_dna(filename, outfilename):
    with open(filename) as fasta_file:
        identifiers = []
        sequences= []
        for seq_record in SeqIO.parse(fasta_file, 'fasta'):  # (generator)
            identifiers.append(seq_record.id)
            sequences.append(seq_record.seq)
    df=np.column_stack((identifiers,sequences))
    df=pd.DataFrame(data=df,columns=['id','seq'])
    a=lambda x: str(x).split(',')
    b = lambda x: x[0]
    df.seq = df.seq.apply(a)
    df.seq = df.seq.apply(b)
    print(df.shape)
    df = df[df.seq.str.len() > 300]
    print(df.shape)
    for i in range(df.shape[0]):
        sequence=[]
        for j in range(0, len(df.iloc[i,1]), 3):
            sequence.append(df.iloc[i,1][j:j+3])
        df.iloc[i,1]=sequence
    text_file = open(outfilename, "w")
    for i in range(df.shape[0]):
        sequence=' '.join(df.iloc[i,1])
        #text_file.write(df.iloc[i,0])
        text_file.write(sequence)
        text_file.write('\n')
    text_file.close()

# function to merge cDNA libraries into one library
def merge_dna(filenamelist, outfilename):
    with open(outfilename, 'w') as outfile:
        for fname in filenamelist:
            with open(fname) as infile:
                outfile.write(infile.read())
