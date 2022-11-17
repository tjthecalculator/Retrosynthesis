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
    [atom.ClearProp('molAtomMapNumber') for atom in Mol.GetAtoms() if atom.HasProp('molAtomMapNumber')]
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
    for atom in Mol.GetAtoms():
        if atom.HasProp('molAtomMapNumber'):
           atom.SetIsotope(int(atom.GetProp('molAtomMapNumber')))

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

def ClearIsotope(Mol):
    [atom.SetIsotope(0) for atom in Mol.GetAtoms()]

def GetChangedAtoms(Reactants, Products):
    ProductAtoms, ProductAtomTags = GetTaggedAtomsFromMols(Products)
    print('Product contain {} tagged atoms'.format(len(ProductAtoms)))
    print('Product contain {} unique atom number'.format(len(ProductAtomTags)))
    ReactantAtoms, ReactantAtomTags = GetTaggedAtomsFromMols(Reactants)
    if len(set(ProductAtomTags)) != len(set(ReactantAtomTags)):
        print('WARNING: Different atom tag appear in reactants and product')
    if len(ProductAtoms) != len(ReactantAtoms):
        print('WARNING: Total number of tagged atoms differ, stoichometry != 1?')
    ChangedAtoms    = []
    ChangedAtomTags = []
    for i, producttag in enumerate(ProductAtomTags):
        for j, reactanttag in enumerate(ReactantAtomTags):
            if reactanttag != producttag:
                continue
            if reactanttag not in ChangedAtomTags:
                if AreAtomDifferent(ProductAtoms[i], ReactantAtoms[j]):
                    ChangedAtoms.append(ReactantAtoms[j])
                    ChangedAtomTags.append(reactanttag)
                    break
                if ProductAtomTags.count(reactanttag) > 1:
                    ChangedAtoms.append(ReactantAtoms[j])
                    ChangedAtomTags.append(reactanttag)
                    break
    for j, reactanttag in enumerate(ReactantAtomTags):
        if reactanttag not in ChangedAtomTags:
            if reactanttag not in ProductAtomTags:
                ChangedAtoms.append(ReactantAtoms[j])
                ChangedAtomTags.append(reactanttag)
    TetraAtoms = GetTetrahedralAtoms(Reactants, Products)
    print('Found {} atom-mapped tetrahedral atoms that have chirality specified at least partially'.format(len(TetraAtoms)))
    [SetIsotope2EqualMapNumber(reactant) for reactant in Reactants]
    [SetIsotope2EqualMapNumber(product) for product in Products]
    for (atomtag, ar, ap) in TetraAtoms:
        print('For atom tag {}'.format(atomtag))
        print('Reactant:    {}'.format(ar.GetChiralTag()))
        print('Product:     {}'.format(ap.GetChiralTag()))
        if atomtag in ChangedAtomTags:
            print('-> atoms have changed (by more than just chirality!)')
        else:
            Unchanged = CheckTetrahedralCenterEquivalent(ar, ap) and Chem.rdchem.ChiralType.CHI_UNSPECIFIED not in [ar.GetChiralTag(), ap.GetChiralTag()]
            if Unchanged:
                print('-> atoms confirmed to have same chirality, no change')
            else:
                print('-> atom changed chirality!!')
                TetraAdjust2RXN = False
                for neighbor in ap.GetNeighbors():
                    if neighbor.HasProp('molAtomMapNumber'):
                        if neighbor.GetProp('molAtomMapNumber') in ChangedAtomTags:
                            TetraAdjust2RXN = True
                            break
                if TetraAdjust2RXN:
                    print('-> atom adjust to reaction center, now included.')
                    ChangedAtoms.append(ar)
                    ChangedAtomTags.append(atomtag)
                else:
                    print('-> adj far from reaction center, not including')
    [ClearIsotope(reactant) for reactant in Reactants]
    [ClearIsotope(product) for product in Products]
    print('{} tagged atoms in reactants change 1-atom properties'.format(len(ChangedAtomTags)))
    for smarts in [atom.GetSmarts() for atom in ChangedAtoms]:
        print('{}'.format(smarts))
    return ChangedAtoms, ChangedAtomTags

def GetSpecialGroups(Mol):
    GroupTemplates = [
        (range(3), '[OH0,SH0]=C[O,Cl,I,Br,F]',),
        (range(3), '[OH0,SH0]=CN',),
        (range(4), 'S(O)(O)Cl',),
        (range(3), 'B(O)O',),
        ((0,), '[Si](C)(C)C'),
        ((0,), '[Si](OC)(OC)(OC)'),
        (range(3), '[N;H0;$(N-[#6]);D2]-,=[N;D2]-,=[N;D1]',),
        (range(8), 'O=C1N([Br,I,F,Cl])C(=O)CC1',),
        (range(11), 'Cc1ccc(S(=O)(=O)O)cc1'),
        ((7,), 'CC(C)(C)OC(=O)[N]'),
        ((4,), '[CH3][CH0]([CH3])([CH3])O'),
        (range(2), '[C,N]=[C,N]',),
        (range(2), '[C,N]#[C,N]',),
        ((2,), 'C=C-[*]',),
        ((2,), 'C#C-[*]',),
        ((2,), 'O=C-[*]',),
        ((3,), 'O=C([CH3])-[*]'),
        ((3,), 'O=C([O,N])-[*]',),
        (range(4), 'ClS(Cl)=O',),
        (range(2), '[Mg,Li,Zn,Sn][Br,Cl,I,F]',),
        (range(3), 'S(O)(O)',),
        (range(2), 'N~N',),
        ((1,), '[!#6;R]@[#6;R]',),
        ((2,), '[a!c]:a:a',),
        ((0,), '[B,C](F)(F)F'),
        ]
    GroupTemplates += [
        ((1,2), '[*]/[CH]=[CH]/[*]'),
        ((1,2), '[*]/[CH]=[CH]\[*]'),
        ((1,2), '[*]/[CH]=[CH0]([*])\[*]'),
        ((1,2), '[*]/[D3;H1]=[!D1]'),
        ]
    Groups = []
    for (addIfMatch, template) in GroupTemplates:
        Matches = Mol.GetSubstructMatches(Chem.MolFromSmarts(template), useChirality=True)
        for match in Matches:
            addIf = []
            for patternIdx, atomIdx in enumerate(match):
                if patternIdx in addIfMatch:
                    addIf.append(atomIdx)
                Groups.append((addIf, match))
    return Groups

def ConvertAtom2Wildcard(Atom):
    if Atom.GetDegree() == 1:
        Symbol = '[' + Atom.GetSymbol() + ';D1;H{}'.format(Atom.GetTotalNumHs())
        if Atom.GetFormalCharge() != 0:
            Charges = re.search('([-+]+[1-9]?)', Atom.GetSmarts())
            Symbol  = Symbol.replace(';D1', ';{};D1'.format(Charges.group()))
    else:
        Symbol = '['
        if Atom.GetAtomicNum() != 6:
            Symbol += '#{};'.format(Atom.GetAtomicNum())
            if Atom.GetIsAromatic():
                Symbol += 'a;'
        elif Atom.GetIsAromatic():
            Symbol += 'c;'
        else:
            Symbol += 'C;'
        if Atom.GetFormalCharge() != 0:
            Charges = re.search('([-+]+[1-9]?)', Atom.GetSmarts())
            if Charges:
                Symbol += Charges.group() + ';'
        if Symbol[-1] == ';':
            Symbol = Symbol[-1]
    Label = re.search('\:[0-9]+\]', Atom.GetSmarts())
    if Label:
        Symbol += Label.group()
    else:
        Symbol += ']'
    if Symbol != Atom.GetSmarts():
        print('Improved generality of atom SMARTS {} -> {}'.format(Atom.GetSmarts(), Symbol))
    return Symbol

def ExpandAtoms2UseAtom(Mol, Atom2Use, AtomIdx, Groups=[], SymbolReplacements=[]):
    FoundInGroup = False
    for group in Groups:
        if int(AtomIdx) in group[0]:
            print('Adding group due to match')
            try:
                print('Match from molAtomMapNum {}'.format(Mol.GetAtomWithIdx(AtomIdx).GetProp('molAtomMapNumber')))
            except KeyError:
                pass
            for idx in group[1]:
                if idx not in Atom2Use:
                    Atom2Use.append(idx)
                    SymbolReplacements.append((idx, ConvertAtom2Wildcard(Mol.GetAtomWithIdx(idx))))
            FoundInGroup = True
    if FoundInGroup:
        return Atom2Use, SymbolReplacements
    if AtomIdx in Atom2Use:
        return Atom2Use, SymbolReplacements
    Atom2Use.append(AtomIdx)
    SymbolReplacements.append((AtomIdx, ConvertAtom2Wildcard(Mol.GetAtomWithIdx(AtomIdx))))
    return Atom2Use, SymbolReplacements

def ExpandAtoms2Use(Mol, Atom2Use, Groups=[], SymbolReplacements=[]):
    NewAtom2Use = Atom2Use[:]
    for atom in Mol.GetAtoms():
        if atom.GetIdx() not in Atom2Use:
            continue
        for group in Groups:
            if int(atom.GetIdx()) in group[0]:
                print('Adding group due to match')
                try:
                    print('Match from molAtomMapNum {}'.format(atom.GetProp('molAtomMapNumber')))
                except KeyError:
                    pass
                for idx in group[1]:
                    if idx not in Atom2Use:
                        NewAtom2Use.append(idx)
                        SymbolReplacements.append((idx, ConvertAtom2Wildcard(Mol.GetAtomWithIdx(idx))))
        for neighbor in atom.GetNeighbors():
            NewAtom2Use, SymbolReplacements = ExpandAtoms2UseAtom(Mol, NewAtom2Use, neighbor.GetIdx(), Groups, SymbolReplacements)
    return NewAtom2Use, SymbolReplacements

def ReassignAtomMapping(Transform):
    AllLabel        = re.findall('\:([0-9]+)\]', Transform)
    Replacements    = []
    ReplacementDict = {}
    Counter         = 1
    for label in AllLabel:
        if label not in ReplacementDict:
            ReplacementDict[label] = str(Counter)
            Counter += 1
        Replacements.append(ReplacementDict[label])
    TransformNewmaps = re.sub('\:[0-9]+\]', lambda match: (':' + Replacements.pop(0) + ']'), Transform)
    return TransformNewmaps

