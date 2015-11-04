# -*- coding: utf-8 -*-
"""
Created on Tue Nov 03 17:05:22 2015

@author: Chris
"""

from rdkit import Chem
from rdkit.Chem import ChemicalFeatures

from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import Draw
from rdkit.Chem.Draw import SimilarityMaps

from rdkit.Chem.Pharm2D import Generate
from rdkit.Chem import BRICS
supp = Chem.SmilesMolSupplier('data/srpk1-diverse.smi')
import search2d


from time import time
import sys


fdefName = 'data/minimalFeatures.fdef'
featFactory = ChemicalFeatures.BuildFeatureFactory(fdefName)


from rdkit.Chem.Pharm2D.SigFactory import SigFactory
sigFactory = SigFactory(featFactory,minPointCount=2,maxPointCount=3, trianglePruneBins=False)
sigFactory.SetBins([(0,2),(2,5),(5,8)])
sigFactory.Init()
sigFactory.GetSigSize()


def similarityMeasure(fps, neg, mol2):
        
    
    fps2 = Generate.Gen2DFingerprint(mol2,sigFactory)
    
    similarityPos = DataStructs.FingerprintSimilarity(fps,fps2, metric=DataStructs.TanimotoSimilarity)
    similarityNeg = DataStructs.FingerprintSimilarity(neg,fps2, metric=DataStructs.TanimotoSimilarity)
#    if similarityPos>=0.75:
            
    print Chem.MolToSmiles(mol2), similarityPos,  similarityNeg
    return similarityPos , similarityPos-similarityNeg   



if __name__ == "__main__":



    allfrags=set()
    for m in supp:
        pieces = BRICS.BRICSDecompose(m)
        allfrags.update(pieces)
        
    print len(allfrags)

    currtime = time()
    #make new molecules from fragments
    import random
    random.seed(127)
    fragms = [Chem.MolFromSmiles(x) for x in allfrags]
    ms = BRICS.BRICSBuild(fragms)

    prods = [ms.next() for x in range(10000)]
    #clean up generated molecules
    for prod in prods:
        prod.UpdatePropertyCache(strict=False)
        
    #srpin340 is a low affinity but selective SRPK1 inhibitor        
        
    srpin340Mol = Chem.MolFromSmiles('C1CCN(CC1)C2=C(C=C(C=C2)C(F)(F)F)NC(=O)C3=CC=NC=C3')
    srpin340fps = Generate.Gen2DFingerprint(srpin340Mol,sigFactory)

    #sphinx is a higher affinity but selective SRPK1 inhibitor 
    
    sphinxMol = Chem.MolFromSmiles('C1(=CC=C(C(=C1)N(C(=O)C2=CC=C(O2)C)[H])N3CCOCC3)C(F)(F)F')
    sphinxfps = Generate.Gen2DFingerprint(sphinxMol,sigFactory)
    
    #CHEMBL2000345 is a broad spectrum kinase inhibitor but also has highest 
    #affinity to srpk1 of all ligands
    #not selective - try to do negative design
    
    
    CHEMBL2000345Mol = Chem.MolFromSmiles('COc1cc(C=C(C#N)c2nc3cc(C)ccc3[nH]2)c(Br)cc1O')
    CHEMBL2000345fps = Generate.Gen2DFingerprint(CHEMBL2000345Mol,sigFactory)    

    # assess new compounds by pharmacophore similarity

    pos_x=[]
    neg_y=[]
    for mol2 in prods:
        pos, neg = similarityMeasure(sphinxfps, CHEMBL2000345fps, mol2)
        if pos>=0.0:
            pos_x.append(pos)
            neg_y.append(neg)
        sys.stdout.flush()

    import pylab as pl
 
    fig,ax = pl.subplots()
    ax.scatter(pos_x, neg_y, alpha=0.3)
    ax.set_xlabel('Positive')
    ax.set_ylabel('Negative')

    ax.legend(loc="lower right")
    
    fig.show()


