# -*- coding: utf-8 -*-
"""
Created on Mon Nov 02 20:49:23 2015

@author: Chris
"""

from rdkit import Chem
from rdkit.Chem import ChemicalFeatures

from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import Draw
from rdkit.Chem.Draw import SimilarityMaps

from rdkit.Chem.Pharm2D import Generate

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
    if similarityPos>=0.75:
            
        print mol2.GetProp('_Name'), Chem.MolToSmiles(mol2), similarityPos,  similarityNeg
    return similarityPos , similarityPos-similarityNeg   
        
        

from multiprocessing import Pool
from time import time
import sys

if __name__ == "__main__":

    currtime = time()
    
    suppl = Chem.SmilesMolSupplier('data/11_p0-test773.smi')
    po = Pool(4)
    
    srpin340Mol = Chem.MolFromSmiles('C1CCN(CC1)C2=C(C=C(C=C2)C(F)(F)F)NC(=O)C3=CC=NC=C3')
    srpin340fps = Generate.Gen2DFingerprint(srpin340Mol,sigFactory)
    
    #CHEMBL2000345 is a broad spectrum kinase inhibitor but also has highest 
    #affinity to srpk1 of all ligands
    #not selective - try to do negative design
    
    
    CHEMBL2000345Mol = Chem.MolFromSmiles('COc1cc(C=C(C#N)c2nc3cc(C)ccc3[nH]2)c(Br)cc1O')
    CHEMBL2000345fps = Generate.Gen2DFingerprint(CHEMBL2000345Mol,sigFactory)    
    pos_x=[]
    neg_y=[]
    for mol2 in suppl:
        pos, neg = similarityMeasure(srpin340fps, CHEMBL2000345fps, mol2)
        if pos>=0.6:
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

    
#    res = po.map_async(similarityMeasure, (srpin340fps, CHEMBL2000345fps, ((mol2) for mol2 in suppl)))
#    print res.get()
#    
    print '2: parallel: time elapsed:', time() - currtime
    
    
    