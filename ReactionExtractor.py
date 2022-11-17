import re
from copy import deepcopy
from rdkit import Chem
from rdkit.Chem import AllChem

USE_STEREOCHEMISTRY                   = True
MAXIMUM_NUMBER_UNMAPPED_PRODUCT_ATOMS = 5
INCLUDE_ALL_UNMAPPED_REACTANT_ATOMS   = True

def MolsFromSMILESList(List):
    Mols = []
    for smiles in List:
        if not smiles:
            continue
        Mols.append(Chem.MolFromSmiles(smiles))
    return Mols

def ReplaceDeuterated(SMILES):
    return re.sub('\[2H]\]', r'[H]', SMILES)

def ClearMapNumber(Mol):
    [a.ClearProp('molAtomMapNumber') for a in Mol.GetAtoms() if a.HasProp('molAtomMapNumber')]
    return Mol

def GetTaggedAtomsFromMol(Mol):
    Atoms    = []
    AtomTags = []
    for atom in Mol.GetAtoms():
        if atom.HasProp('molAtomMapNumber'):
            Atoms.append(atom)
            AtomTags.append(str(atom.GetProp('molAtomMapNumber')))
    return Atoms, AtomTags

def GetTaggedAtomsFromMols(Mols):
    Atoms    = []
    AtomTags = []
    for mol in Mols:
        NewAtoms, NewAtomTags = GetTaggedAtomsFromMol(mol)
        Atoms    += NewAtoms
        AtomTags += NewAtomTags
    return Atoms, AtomTags

def Bond2Label(Bond):
    Atom1Label = str(Bond.GetBeginAtom().GetAtomicNum())
    Atom2Label = str(Bond.GetEndAtom().GetAtomicNum())
    if Bond.GetBeginAtom().HasProp('molAtomMapNumber'):
        Atom1Label += Bond.GetBeginAtom().GetProp('molAtomMapNumber')
    if Bond.GetEndAtom().HasProp('molAtomMapNumber'):
        Atom2Label += Bond.GetEndAtom().GetProp('molAtomMapNumber')
    Atoms = sorted([Atom1Label, Atom2Label])
    return '{}{}{}'.format(Atoms[0], Bond.GetSmarts(), Atoms[1])

def AreAtomDifferent(Atom1, Atom2):
    if Atom1.GetAtomicNum() != Atom2.GetAtomicNum():
        return True
    if Atom1.GetTotalNumHs() != Atom2.GetTotalNumHs():
        return True
    if Atom1.GetFormalCharge() != Atom2.GetFormalCharge():
        return True
    if Atom1.GetNumRadicalElectrons() != Atom2.GetNumRadicalElectrons():
        return True
    if Atom1.GetIsAromatic() != Atom2.GetIsAromatic():
        return True
    Bonds1 = sorted([Bond2Label(bond) for bond in Atom1.GetBonds()])
    Bonds2 = sorted([Bond2Label(bond) for bond in Atom2.GetBonds()])
    if Bonds1 != Bonds2:
        return True

    return False

def FindMapNumber(Mol, MapNumber):
    return [(a.GetIdx(), a) for a in Mol.GetAtoms() if a.HasProp('molAtomMapNumber') and a.GetProp('molAtomMapNumber') == str(MapNumber)][0]

def GetTetrahedralAtoms(Reactants, Products):
    TetrahedralAtoms = []
    for reactant in Reactants:
        for ar in reactant.GetAtoms():
            if not ar.HasProp('molAtomMapNumber'):
                continue
            AtomTag = ar.GetProp('molAtomMapNumber')
            ir      = ar.GetIdx()
            for product in Products:
                try:
                    (ip, ap) = FindMapNumber(product, AtomTag)
                    if ar.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED or ap.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                        TetrahedralAtoms.append((AtomTag, ar, ap))
                except IndexError:
                    pass

    return TetrahedralAtoms

def SetIsotope2EqualMapNumber(Mol):
    for a in Mol.GetAtoms():
        if a.HasProp('molAtomMapNumber'):
           a.SetIsotope(int(a.GetProp('molAtomMapNumber')))

def GetFragmentAroundTetrahedralCenter(Mol, Idx):
    Idx2Include = [Idx]
    for neighbor in Mol.GetAtomWithIdx(Idx).GetNeighbors():
        Idx2Include.append(neighbor.GetIdx())
    Symbols = ['[{}{}]'.format(atom.GetIsotope(), atom.GetSymbol()) if atom.GetIsotope() != 0 else '[#{}]'.format(atom.GetAtomicNum()) for atom in Mol.GetAtoms()]
    return Chem.MolFragmentToSmiles(Mol, Idx2Include, isomericSmiles=True, atomSymbols=Symbols, allBondsExplicit=True, allHsExplicit=True)

def CheckTetrahedralCenterEquivalent(Atom1, Atom2):
    Atom1Fragments    = GetFragmentAroundTetrahedralCenter(Atom1.GetOwningMol(), Atom1.GetIdx())
    Atom1Neighborhood = Chem.MolFromSmiles(Atom1Fragments, sanitize=False)
    for matchedIds in Atom2.GetOwningMol().GetSubstructMatches(Atom1Neighborhood, useChirality=True):
        if Atom2.GetIdx() in matchedIds:
            return True
    return False
