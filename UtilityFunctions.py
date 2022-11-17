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
        Outcome   = RDKit Products (tuple)
    """
    try:
        Outcome = Reaction.RunReactants(Reactants)
    except ValueError:
        print("Reaction Error!!!! Something wrong in reactant or reaction")
        Outcome = []

    if Outcome:
        for outcome in Outcome:
            for mol in outcome:
                try:
                    Chem.SanitizeMol(mol)
                except ValueError:
                    print("Problem with sanitizing molecule, passing original SMILES format")
                    continue

    return Outcome

def RewardFunction(Mols):
    """
    Received RDKit molecule then calculating reward score based on is there any similar molecule in dataset.

    Input:
        Mols = RDKit molecules (tuple)

    Output:
        Reward = Reward score
    """
    Reward = 0
    for mol in Mols:
        pass

def ActionValue(Node, Exploration=3.0):
    """
    
    """
    Q = Node.GetActionValue()/Node.GetNumVisited()
    U = Exploration*Node.GetPriorProbability()*np.sqrt(Node.GetParentNode().GetNumVisited())/(1 + Node.GetNumVisited())
    return Q + U

def UpdateActionValue(ParentNode, LMAX=10, DAMPING=0.99):
    pass