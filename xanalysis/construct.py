from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import pandas as pd
import numpy as np
import random

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

# make protein subsample with a- and b-mers as words
def kmerize_prot(filename, outfilename, n=1000, a=1, b=2):
    i=0
    with open(outfilename,'w') as outfile:
        while i<n:
            outline=[]
            line = random.choice(list(open(filename)))
            line = line.replace('\n','')
            if len(line)%2==1:
                line=line[:-1]
            j=0
            while j<len(line):
                choice = random.choice([a,b])
                if choice==a:
                    outline.append(line[j:j+a])
                    j+=a
                elif choice==b:
                    outline.append(line[j:j+b])
                    j+=b
            outline=str(' '.join(outline))
            i+=1
            outfile.write(outline+'\n')

# construct feature matrix from amino acid sequence
def featurize_prot(filename,idColumn,seqColumn,label):
    df=pd.read_csv(filename)
    df[seqColumn]=df[seqColumn].apply(lambda x: x.replace('\r\n',''))
    idx=df[idColumn]
    seqs=df[seqColumn].values.tolist()
    weight = []
    isoelectric = []
    helix =[]
    turn=[]
    sheet=[]
    aromaticity=[]
    A=[]
    C=[]
    D=[]
    E=[]
    F=[]
    G=[]
    H=[]
    I=[]
    K=[]
    L=[]
    M=[]
    N=[]
    P=[]
    Q=[]
    R=[]
    S=[]
    T=[]
    V=[]
    W=[]
    Y=[]
    flex1=[]
    flex2=[]
    flex3=[]
    ext=[]
    lab=[]
    dic2=['A','C',
               'D','E','F','G','H','I','K','L','M','N','P',
               'Q','R','S','T','V','W','Y']
    final_dic={}
    for seq in seqs:
        seq=str(seq)
        analyse=ProteinAnalysis(seq)
        weight.append(analyse.molecular_weight())
        isoelectric.append(analyse.isoelectric_point())
        h, t, s = analyse.secondary_structure_fraction()
        helix.append(h)
        turn.append(t)
        sheet.append(h)
        aromaticity.append(analyse.aromaticity())
        dic = analyse.get_amino_acids_percent()
        A.append(dic['A'])
        C.append(dic['C'])
        D.append(dic['D'])
        E.append(dic['E'])
        F.append(dic['F'])
        G.append(dic['G'])
        H.append(dic['H'])
        I.append(dic['I'])
        K.append(dic['K'])
        L.append(dic['L'])
        M.append(dic['M'])
        N.append(dic['N'])
        P.append(dic['P'])
        Q.append(dic['Q'])
        R.append(dic['R'])
        S.append(dic['S'])
        T.append(dic['T'])
        V.append(dic['V'])
        W.append(dic['W'])
        Y.append(dic['Y'])
        flexy=analyse.flexibility()
        _s =round((len(flexy)/3))+1
        flex1.append(np.mean(flexy[:_s]))
        flex2.append(np.mean(flexy[_s:2*_s]))
        flex3.append(np.mean(flexy[2*_s:]))
        ext.append(analyse.molar_extinction_coefficient()[0])
        lab.append(label)
        for o in dic2:
            for j in dic2:
                final_dic.setdefault(o+j,[])
                final_dic[o+j].append(seq.count(o+j)/(len(seq)/2))
    df3 = pd.DataFrame.from_dict(final_dic)
    df2 = np.column_stack((weight,isoelectric,helix,
                          turn,sheet,aromaticity,A,C,D,E,F,G,H,
                          I,K,L,M,N,P,Q,R,S,T,V,W,Y,flex1,flex2,
                           flex3,ext,lab))
    columns = ['MolWeight','IsoelectricPoint',
               'Helix%','Turn%','Sheet%','Aromaticity','A%','C%',
               'D%','E%','F%','G%','H%','I%','K%','L%','M%','N%','P%',
               'Q%','R%','S%','T%','V%','W%','Y%','MeanFlexibility1',
               'MeanFlexibility2','MeanFlexibility3',
               'ExtinctionCoef','FP']
    df2 = pd.DataFrame(data=df2,columns=columns)
    df4 = pd.concat([df2,df3],axis=1)
    return df4
