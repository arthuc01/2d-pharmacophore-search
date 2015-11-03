# -*- coding: utf-8 -*-
"""
Created on Mon Nov 02 20:49:23 2015

@author: Chris
"""

from rdkit import Chem
from rdkit.Chem import ChemicalFeatures

from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols

fdefName = 'data/minimalFeatures.fdef'
featFactory = ChemicalFeatures.BuildFeatureFactory(fdefName)


from rdkit.Chem.Pharm2D.SigFactory import SigFactory
sigFactory = SigFactory(featFactory,minPointCount=2,maxPointCount=3, trianglePruneBins=False)
sigFactory.SetBins([(0,2),(2,5),(5,8)])
sigFactory.Init()
sigFactory.GetSigSize()

from rdkit.Chem.Pharm2D import Generate
mol = Chem.MolFromSmiles('C1CCN(CC1)C2=C(C=C(C=C2)C(F)(F)F)NC(=O)C3=CC=NC=C3')

fps = Generate.Gen2DFingerprint(mol,sigFactory)

from rdkit.Chem import Draw
from rdkit.Chem.Draw import SimilarityMaps

suppl = Chem.SmilesMolSupplier('data/11_p0.smi')

for mol2 in suppl:
        
    
    fps2 = Generate.Gen2DFingerprint(mol2,sigFactory)
    
    similarity = DataStructs.FingerprintSimilarity(fps,fps2, metric=DataStructs.TanimotoSimilarity)
    
    if similarity>=0.6:
        print similarity, mol2.GetProp('_Name'), Chem.MolToSmiles(mol2)