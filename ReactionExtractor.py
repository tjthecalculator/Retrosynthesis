import re
from numpy.random import shuffle
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

def GetStrictSMARTSAtom(Atom):
    Symbol = Atom.GetSmarts()
    if Atom.GetSymbol() == 'H':
        Symbol = '[#1]'
    if '[' not in Symbol:
        Symbol = '[' + Symbol + ']'
    if USE_STEREOCHEMISTRY:
        if Atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
            if '@' not in Symbol:
                if Atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
                    tag = '@'
                elif Atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
                    tag = '@@'
                if ':' in Symbol:
                    Symbol = Symbol.replace(':', ';{}:'.format(tag))
                else:
                    Symbol = Symbol.replace(']', ';{}]'.format(tag))
    if 'H' not in Symbol:
        H_Symbol = 'H{}'.format(Atom.GetTotalNumHs())
        if ':' in Symbol:
            Symbol = Symbol.replace(':', ';{}:'.format(H_Symbol))
        else:
            Symbol = Symbol.replace(']', ';{}]'.format(H_Symbol))
    if ':' in Symbol:
        Symbol = Symbol.replace(':', ';D{}:'.format(Atom.GetDegree()))
    else:
        Symbol = Symbol.replace(']', ';D{}]'.format(Atom.GetDegree()))
    if '+' not in Symbol and '-' not in Symbol:
        Charge       = Atom.GetFormalCharge()
        ChargeSymbol = '+' if (Charge >= 0) else '-'
        ChargeSymbol += '{}'.format(abs(Charge))
        if ':' in Symbol:
            Symbol = Symbol.replace(':', ';{}:'.format(ChargeSymbol))
        else:
            Symbol = Symbol.replace(']', ';{}]'.format(ChargeSymbol))
    return Symbol

def ExpandChangedAtomTags(ChangedAtomTags, ReactantFragments):
    Expansion = []
    AtomTagsReactantFragments = re.findall('\:([0-9]+)\]', ReactantFragments)
    for atomTag in AtomTagsReactantFragments:
        if atomTag not in ChangedAtomTags:
            Expansion.append(atomTag)
    print('After building reactant fragment, additional labels included: {}'.format(Expansion))
    return Expansion

def GetFragmentforChangedAtoms(Mols, ChangedAtomTags, Radius=0, Category='Reactants', Expansion=[]):
    Fragments   = ''
    MolsChanged = []
    for mol in Mols:
        SymbolReplacements = []
        if Category == 'Reactants':
            Groups = GetSpecialGroups(mol)
        else:
            Groups = []
        Atom2Use = []
        for atom in mol.GetAtoms():
            if ':' in atom.GetSmarts():
                if atom.GetSmarts().split(':')[1][:-1] in ChangedAtomTags:
                    Atom2Use.append(atom.GetIdx())
                    Symbol = GetStrictSMARTSAtom(atom)
                    if Symbol != atom.GetSmarts():
                        SymbolReplacements.append((atom.GetIdx(), Symbol))
                    continue
        if INCLUDE_ALL_UNMAPPED_REACTANT_ATOMS and len(Atom2Use) > 0:
            if Category == 'Reactants':
                for atom in mol.GetAtoms():
                    if not atom.HasProp('molAtomMapNumber'):
                        Atom2Use.append(atom.GetIdx())
        for k in range(Radius):
            Atom2Use, SymbolReplacements = ExpandAtoms2Use(mol, Atom2Use, Groups, SymbolReplacements)
        if Category == 'Products':
            if Expansion:
                for atom in mol.GetAtoms():
                    if ':' not in atom.GetSmarts():
                        continue
                    Label = atom.GetSmarts().split(':')[1][:-1]
                    if Label in Expansion and Label not in ChangedAtomTags:
                        Atom2Use.append(atom.GetIdx())
                        SymbolReplacements.append((atom.GetIdx(), ConvertAtom2Wildcard(atom)))
                        print('Expanded label {} to wildcard in products'.format(Label))
            for atom in mol.GetAtoms():
                if not atom.HasProp('molAtomMapNumber'):
                    Atom2Use.append(atom.GetIdx())
                    Symbol = GetStrictSMARTSAtom(atom)
                    SymbolReplacements.append((atom.GetIdx(), Symbol))
        Symbols = [atom.GetSmarts() for atom in mol.GetAtoms()]
        for (i, symbol) in SymbolReplacements:
            Symbols[i] = symbol
        if not Atom2Use:
            continue
        TetraConsistent = False
        NumTetraFlips   = 0
        while not TetraConsistent and NumTetraFlips < 100:
            MolCopy = deepcopy(mol)
            [atom.ClearProp('molAtomMapNumber') for atom in MolCopy.GetAtoms()]
            ThisFragments   = AllChem.MolFragmentToSmiles(MolCopy, Atom2Use, atomSymbols=Symbols, allHsExplicit=True, isomericSmiles=USE_STEREOCHEMISTRY, allBondExplicit=True)
            ThisFragmentMol = AllChem.MolFromSmarts(ThisFragments)
            TetraMapNums    = []
            for atom in ThisFragmentMol.GetAtoms():
                if atom.HasProp('molAtomMapNumber'):
                    atom.SetIsotope(int(atom.GetProp('molAtomMapNumber')))
                    if atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                        TetraMapNums.append(atom.GetProp('molAtomMapNumber'))
            Map2ID = {}
            for atom in mol.GetAtoms():
                if atom.HasProp('molAtomMapNumber'):
                    atom.SetIsotope(int(atom.GetProp('molAtomMapNumber')))
                    Map2ID[atom.GetProp('molAtomMapNumber')] = atom.GetIdx()
            TetraConsistent = True
            AllMatchedIds   = []
            FragmentSMILES  = Chem.MolToSmiles(ThisFragmentMol)
            if FragmentSMILES.count('.') > 5:
                break
            for matchedIds in mol.GetSubstructMatches(ThisFragmentMol, useChirality=True):
                AllMatchedIds.extend(matchedIds)
            shuffle(TetraMapNums)
            for tetraMapNum in TetraMapNums:
                print('Checking consistency of tetrahedral {}'.format(tetraMapNum))
                if Map2ID[tetraMapNum] not in AllMatchedIds:
                    TetraConsistent = False
                    print('@@@@@@@@@@@@@@@@ FRAGMENT DOES NOT MATCH PARENT MOL @@@@@@@@@@@@@@')
                    print('@@@@@@@@@@@@@@@@ FILPPING CHIRALITY SYMBOL NOW      @@@@@@@@@@@@@@')
                    PreviousSymbol = Symbols[Map2ID[tetraMapNum]]
                    if '@@' in PreviousSymbol:
                        Symbol = PreviousSymbol.replace('@@', '@')
                    elif '@' in PreviousSymbol:
                        Symbol = PreviousSymbol.replace('@', '@@')
                    else:
                        raise ValueError('Need to modify symbol of tetra atom without @ or @@???')
                    Symbols[Map2ID[tetraMapNum]] = Symbol
                    NumTetraFlips += 1
                    break
            for atom in mol.GetAtoms():
                atom.SetIsotope(0)
        if not TetraConsistent:
            raise ValueError('Could not find consistent tetrahedral mapping, {} center'.format(len(TetraMapNums)))
        [atom.ClearProp('molAtomMapNumber') for atom in mol.GetAtoms()]
        ThisFragment = AllChem.MolFragmentToSmiles(mol, Atom2Use, atomSymbols=Symbols, allHsExplicit=True, isomericSmiles=USE_STEREOCHEMISTRY, allBondExplicit=True)
        Fragments += '(' + ThisFragment + ').'
        MolsChanged.append(Chem.MolToSmiles(ClearMapNumber(Chem.MolFromSmiles(Chem.MolToSmiles(mol, True))), True))
    IntraOnly = (1 == len(MolsChanged))
    DimerOnly = (1 == len(set(MolsChanged))) and (len(MolsChanged) == 2)
    return Fragments[:-1], IntraOnly, DimerOnly

def CanonicalTemplate(Template):
    TemplateNolabels     = re.sub('\:[0-9]+\]', Template)
    TemplateNolabelsMols = TemplateNolabels[1:-1].split(').(')
    TemplateMols         = Template[1:-1].split(').(')
    for i in range(len(TemplateMols)):
        NolabelMolFragments  = TemplateNolabelsMols[i].split('.')
        MolFragments         = TemplateMols[i].split('.')
        Sortorder            = [j[0] for j in sorted(enumerate(NolabelMolFragments), key=lambda x:x[1])]
        TemplateNolabelsMols = '.'.join([NolabelMolFragments[j] for j in Sortorder])
        TemplateMols         = '.'.join([MolFragments[j] for j in Sortorder])
    Sortorder = [j[0] for j in sorted(enumerate(TemplateNolabelsMols), key=lambda x:x[1])]
    Template  = '(' + ').('.join([TemplateMols[i] for i in Sortorder]) + ')'
    return Template

def CanonicalTransform(Transform):
    TransformReordered = '>>'.join([CanonicalTemplate(x) for x in Transform.split('>>')])
    return ReassignAtomMapping(TransformReordered)

def GetReactingMols(Mols, ChangedAtomTags):
    MolsSMILES = []
    for mol in Mols:
        for atom in mol.GetAtoms():
            if ':' in atom.GetSmarts():
                if atom.GetSmarts().split(':')[1][:-1] in ChangedAtomTags:
                    SMILES = Chem.MolToSmiles(mol)
                    if SMILES not in MolsSMILES:
                        MolsSMILES.append(SMILES)
                    continue
    return MolsSMILES

def Extractor(Reaction):
    Reactants = MolsFromSMILESList(ReplaceDeuterated(Reaction['Reactants']))
    Products  = MolsFromSMILESList(ReplaceDeuterated(Reaction['Products']))
    if None in Reactants or None in Products:
        return {'Reaction_ID':Reaction['ID']}
    try:
        for i in range(len(Reactants)):
            Reactants[i] = AllChem.RemoveHs(Reactants[i])
        for i in range(len(Products)):
            Products[i]  = AllChem.RemoveHs(Products[i])
        [Chem.SanitizeMol(mol) for mol in Reactants + Products]
        [mol.UpdatePropertyCache() for mol in Reactants + Products]
    except Exception as e:
        print(e)
        print('Could not load SMILES or Sanitize')
        print('ID: {}'.format(Reaction['ID']))
        return {'Reaction_ID':Reaction['ID']}
    AreUnmappedProductAtoms = False
    ExtraReactantFragment   = ''
    for product in Products:
        productAtoms = product.GetAtoms()
        if sum([atom.HasProp('molAtomMapNumber') for atom in productAtoms]) < len(productAtoms):
            print('Not all product atoms have atom mapping')
            print('ID: {}'.format(Reaction['ID']))
            AreUnmappedProductAtoms = True
    if AreUnmappedProductAtoms:
        for product in Products:
            productAtoms = product.GetAtoms()
            UnmappedIds  = [atom.GetIdx() for atom in productAtoms if not atom.HasProp('molAtomMapNumber')]
            if len(UnmappedIds) > MAXIMUM_NUMBER_UNMAPPED_PRODUCT_ATOMS:
                print('Skip this example - too many unmapped product atoms!')
                print('ID: {}'.format(Reaction['ID']))
                return {'Reaction_ID':Reaction['ID']}
            AtomSymbols = ['[{}]'.format(atom.GetSymbol()) for atom in productAtoms]
            BondSymbols = ['~' for bond in product.GetBonds()]
            if UnmappedIds:
                ExtraReactantFragment += AllChem.MolFragmentToSmiles(product, UnmappedIds, allHsExplicit=True, isomericSmiles=USE_STEREOCHEMISTRY, atomSymbols=AtomSymbols, bondSymbols=BondSymbols) + '.'
        if ExtraReactantFragment:
            ExtraReactantFragment = ExtraReactantFragment[:-1]
            print('Extra reactant fragement: {}'.format(ExtraReactantFragment))
        ExtraReactantFragment = '.'.join(sorted(list(set(ExtraReactantFragment.split('.')))))
    if None in Reactants + Products:
        print('Could not parse all molecules in reaction, skipping')
        print('ID: {}'.format(Reaction['ID']))
        return {'Reaction_ID':Reaction['ID']}
    ChangedAtoms, ChangedAtomTags = GetChangedAtoms(Reactants, Products)
    if not ChangedAtomTags:
        print('No atom changed?')
        print('ID: {}'.format(Reaction['ID']))
        return {'Reaction_ID':Reaction['ID']}
    MolsReactants = GetReactingMols(Reactants, ChangedAtomTags)
    MolsProducts  = GetReactingMols(Products, ChangedAtomTags)
    try:
        ReactantFragments, IntraOnly, DimerOnly = GetFragmentforChangedAtoms(Reactants, ChangedAtomTags, Radius=1, Expansion=[], Category='Reactants')
        ProductFragments, _, _ = GetFragmentforChangedAtoms(Products, ChangedAtomTags, Radius=0, Expansion=ExpandChangedAtomTags(ChangedAtomTags, ReactantFragments), Category='Products')
    except ValueError as e:
        print(e)
        print('ID: {}'.format(Reaction['ID']))
        return {'Reaction_ID':Reaction['ID']}
    RXNString         = '{}>>{}'.format(ReactantFragments, ProductFragments)
    RXNCanonical      = CanonicalTransform(RXNString)
    RXNCanonicalSplit = RXNCanonical.split('>>')
    RXNCanonical      = RXNCanonicalSplit[0][1:-1].replace(').(', '.') + RXNCanonicalSplit[1][1:-1].replace(').(', '.')
    ReactantsString   = RXNCanonical.split('>>')[0]
    ProductsString    = RXNCanonical.split('>>')[1]
    RetroCanonical    = ProductsString + '>>' + ReactantsString
    try:
        RXN = AllChem.ReactionFromSmarts(RetroCanonical)
        if RXN.Validate()[1] != 0:
            print('Could not validate reaction successfully')
            print('ID: {}'.format(Reaction['ID']))
            print('Retro Canonical: {}'.format(RetroCanonical))
            return {'Reaction_ID':Reaction['ID']}
    except Exception as e:
        print(e)
        print('ID: {}'.format(Reaction['ID']))
        return {'Reaction_ID':Reaction['ID']}
    Templates = {
        'Products':MolsProducts,
        'Reactants':MolsReactants,
        'Reaction_SMARTS':RetroCanonical,
        'Intra_Only':IntraOnly,
        'Dimer_Only':DimerOnly,
        'Reaction_ID':Reaction['ID'],
        'Necessary_Reagent':ExtraReactantFragment,
        'Spectators':Reaction['Spectators']
        }
    return Templates