from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np

def Fingerprint2Array(Mol, Size):
    """
    Received RDKit single molecule then calculate fingerprint and convert into numpy array according to size.

    Input:
        Mol  = RDKit molecule
        Size = Fingerprint size 

    Output:
        NFP  = Numpy array of fingerprint
    """
    FP = AllChem.GetMorganFingerprintAsBitVect(Mol, radius=2, useChirality=True, nBits=Size)
    NFP = np.zeros((0,))
    Chem.DataStructs.ConvertToNumpyArray(FP, NFP)
    return NFP

def RunReaction(Reaction, Reactants):
    """
    Received RDKit Reaction and Reactant then generate product according to reaction then sanitize the molecule.

    Input:
        Reaction  = RDKit Reaction
        Reactants = RDKit Reactant (tuple)

    Output:
        
    """
    try:
        Outcome = Reaction.RunReactants(Reactants)
    except ValueError:
        print("Reaction Error!!!! Something wrong in reactant or reaction")
        Outcome = []

    if Outcome:
