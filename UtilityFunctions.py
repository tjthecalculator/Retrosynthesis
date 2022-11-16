from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

def Fingerprint2Array(Mol, Size):
    """
    Received RDKit single molecule then calculate fingerprint and convert into numpy array according to size.

    Input:
        Mol  = RDKit molecule
        Size = Fingerprint size 
    """
    FP = AllChem.GetMorganFingerprintAsBitVect(Mol, radius=2, useChirality=True, nBits=Size)
    NFP = np.zeros((0,))
    Chem.DataStructs.ConvertToNumpyArray(FP, NFP)
    return NFP