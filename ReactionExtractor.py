import re
from copy import deepcopy
from rdkit import Chem
from rdkit.Chem import AllChem

USE_STEREOCHEMISTRY                   = True
MAXIMUM_NUMBER_UNMAPPED_PRODUCT_ATOMS = 5
INCLUDE_ALL_UNMAPPED_REACTANT_ATOMS   = True

def Mols_From_SMILESList(List):
    Mols = []
    for smiles in List:
        if not smiles:
            continue
        Mols.append(Chem.MolFromSmiles(smiles))
    return Mols
