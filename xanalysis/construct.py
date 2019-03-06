from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import pandas as pd
import numpy as np
import random


def get_GeneID_list(n=10000,out='list.txt'):
    """ function to generate a list of GeneID
    identifiers to get random proteins from UniProt """
    ls=np.random.choice(2700000,n)
    with open(out,'w') as file:
        for i in ls:
            file.write(str(i)+'\n')
    file.close()


def translate_aa(filename, outputfile):
    """ function to transform cDNA library into protein library """
    with open(filename) as inputfile:
        for line in inputfile:
                line = line.replace(' ','')
                line = line.replace('\n','')
                line = Seq(line, IUPAC.unambiguous_dna)
                line = str(line.translate())
                line = line.replace('*','')
                with open(outputfile,'a') as the_file:
                    the_file.write(line+'\n')


def prep_dna(filename, outfilename):
    """ function to delete fasta ids, short ORFs and
    splits dna sequences into triplets one gene per line """
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

 
def prep_prot(filename,outfilename):
    """ function to delete fasta ids and short proteins, one protein per line """
    fasta_seqs=SeqIO.parse(open(filename),'fasta')
    identifiers=[]
    sequences=[]
    for fasta in fasta_seqs:
        identifiers.append(fasta.id)
        sequences.append(fasta.seq)
    df=np.column_stack((identifiers,sequences))
    df=pd.DataFrame(data=df,columns=['id','seq'])
    a=lambda x: str(x).split(',')
    b = lambda x: x[0]
    df.seq = df.seq.apply(a)
    df.seq = df.seq.apply(b)
    print('Proteins in total:',df.shape)
    df = df[df.seq.str.len() > 100]
    df = df[df.seq.str.len() < 1000]
    print('Proteins between 100-1000 aa:',df.shape)
    text_file=open(outfilename,'w')
    for i in range(df.shape[0]):
        sequence=''.join(df.iloc[i,1])
        text_file.write(sequence)
        text_file.write('\n')
    text_file.close()
    

def merge_dna(filenamelist, outfilename):
    """ function to merge cDNA libraries into one library """
    with open(outfilename, 'w') as outfile:
        for fname in filenamelist:
            with open(fname) as infile:
                outfile.write(infile.read())


def fasta_to_csv(filename,class_label,outfilename):
    """ function to extract sequence from
    fasta files and convert it into a csv file"""
    fasta_seqs=SeqIO.parse(open(filename),'fasta')
    labels=[]
    seqs=[]
    for fasta in fasta_seqs:
        if all(s not in fasta for s in ['B','J','O','U','X','Z']):
            labels.append(int(class_label))
            seqs.append(fasta.seq)
    df=np.column_stack((labels,seqs))
    df=pd.DataFrame(data=df,columns=['Class','Sequence'])
    df.to_csv(outfilename,encoding='utf-8',index=False)


def csv_to_text(filename,seqColumn,outfilename):
    """ function to convert csv file of prot sequences into txt file """
    df=pd.read_csv(filename)
    col_len=df.shape[0]
    with open(outfilename,'w') as file:
        for i in range(col_len):
            line=df.loc[i,seqColumn]
            line=line.strip()
            line=line.replace('\n','').replace('\r','')
            file.write(line+'\n')
        file.close()


def text_to_csv(filename,outfilename):
    """ function to convert text output to csv sequence file """
    seqs=[]
    with open(filename,'r') as file:
        for line in file:
            seqs.append(line.replace(' ','').replace('\n',''))
    df=pd.DataFrame(data=seqs,columns=['Sequence'])
    df.to_csv(outfilename,encoding='utf-8',index=False)


def source_to_dataset(neg_examples,pos_examples,outfilename,
                      neg_seq='Sequence',
                      pos_seq='Sequence'):
    """ takes negative & positive training examples
    and constructs a featurized dataset """
    df_neg = featurize_prot(neg_examples,neg_seq,0)
    df_pos = featurize_prot(pos_examples,pos_seq,1)
    df = pd.concat([df_neg,df_pos],axis=0)
    df.to_csv(outfilename,encoding='utf-8',index=False)


def kmerize_prot(filename, outfilename, n=1000, a=1, b=2,c=3):
    """ make protein subsample with a- and b-mers as words """
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
                choice = random.choice([a,b,c])
                if choice==a:
                    outline.append(line[j:j+a])
                    j+=a
                elif choice==b:
                    outline.append(line[j:j+b])
                    j+=b
                elif choice==c:
                    outline.append(line[j:j+c])
                    j+=c
            outline=str(' '.join(outline))
            i+=1
            outfile.write(outline+'\n')


def featurize_prot(filename,seqColumn,label):
    """ construct feature matrix from amino acid sequence """
    df=pd.read_csv(filename)
    df[seqColumn]=df[seqColumn].apply(lambda x: x.replace('\r\n',''))
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
    instab=[]
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
        instab.append(analyse.instability_index())
        for o in dic2:
            for j in dic2:
                final_dic.setdefault(o+j,[])
                final_dic[o+j].append(seq.count(o+j)/(len(seq)/2))
    df3 = pd.DataFrame.from_dict(final_dic)
    df2 = np.column_stack((weight,isoelectric,helix,
                          turn,sheet,aromaticity,A,C,D,E,F,G,H,
                          I,K,L,M,N,P,Q,R,S,T,V,W,Y,flex1,flex2,
                           flex3,ext,instab,lab))
    columns = ['MolWeight','IsoelectricPoint',
               'Helix%','Turn%','Sheet%','Aromaticity','A%','C%',
               'D%','E%','F%','G%','H%','I%','K%','L%','M%','N%','P%',
               'Q%','R%','S%','T%','V%','W%','Y%','MeanFlexibility1',
               'MeanFlexibility2','MeanFlexibility3',
               'ExtinctionCoef','Instability','Fluor']
    df2 = pd.DataFrame(data=df2,columns=columns)
    df4 = pd.concat([df2,df3],axis=1)
    df4 = df4.loc[df4['Instability']<40.0]
    return df4
