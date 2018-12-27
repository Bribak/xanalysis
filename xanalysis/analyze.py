from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import pandas as pd
import numpy as np

# function to analyze protein properties of a protein library
def analyze(filename):
    instability = []
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
    with open(filename) as inputfile:
        for line in inputfile:
            line = line.replace(' ','')
            line = line.replace('\n','')
            line = Seq(line, IUPAC.unambiguous_dna)
            line = str(line.translate())
            line = line.replace('*','')
            analyse = ProteinAnalysis(line)
            instability.append(analyse.instability_index())
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
    df = np.column_stack((weight,instability,isoelectric,helix,
                          turn,sheet,aromaticity,A,C,D,E,F,G,H,
                          I,K,L,M,N,P,Q,R,S,T,V,W,Y))
    columns = ['MolWeight','Instability','IsoelectricPoint',
               'Helix%','Turn%','Sheet%','Aromaticity','A%','C%',
               'D%','E%','F%','G%','H%','I%','K%','L%','M%','N%','P%',
               'Q%','R%','S%','T%','V%','W%','Y%']
    df = pd.DataFrame(data=df,columns=columns)
    return df
