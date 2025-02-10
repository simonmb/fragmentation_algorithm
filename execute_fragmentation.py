import operator
from rdkit import Chem
from fragmenter import fragmenter
from fragmenter_utils import draw_mol_with_highlights_and_legend
import SMARTS
from pprint import pprint

UNIFAC_SMARTS = SMARTS.UNIFAC.copy()


sorted_group_names_as_in_paper = [
    'Furfural', 'NMP', '(CH3O)2CO', 'MORPH', 'C5H5N', 'CCl4', 'CCl3F', 
    'CClF3', 'CCl2F2', 'DMF', 'C4H4S', 'CHCl3', 'HCCl2F', 'HCClF2', 
    'DMSO', 'DOH', 'ACRY', 'CS2', 'HCOOH', 'CH2Cl2', 'CH3CN', 'H2COCH2',
    'CH3OH', 'CH3NH2', 'CH3SH', 'H2O', 'CCl3', 'CF3', 'CCl2F', 'CClF2',
    'CHCl2', 'HCClF', 'Cl(C=C)', 'CCl2', 'CF2', 'CH2Cl', 'ACCl', 'ACF',
    'CHCl', 'CCl', 'CF', 'I', 'Br', 'CH3OCH2OCO', 'C5H4N', '(CH2O)2CO',
    'C5H3N', 'CON(CH3)2', 'C4H3S', 'OCOCO', '(CH2)2SU', 'CH2CHSU',
    'HCON(CH2)2', 'CONCH3CH2', 'CON(CH2)2', 'C4H2S', 'CH3NO2', 'CH3COO',
    'CONHCH3', 'C2H5O2', 'ACNO2', 'CH2NO2', 'CHNO2', 'CH2COO', 'CONHCH2',
    'C2H4O2', 'NCO', 'HCOO', 'COOH', 'CONH2', 'CH2CN', 'CH3CO', 'H2COCH',
    'COO', 'CH2CO', 'HCOCH', 'COCH', 'CH=O', 'CH3O', 'CH2NH2', 'CH3NH',
    'CH2SH', 'CH3S', 'CH#C', 'CH2=CH', 'SiH2O', 'SiHO', 'SiO', 'ACOH',
    'THF', 'ACNH2', 'CH2O', 'CHO', 'CHNH2', 'CH2NH', 'CHNH', 'CH3N',
    'CH2N', 'CH2S', 'CHS', 'ACCH3', 'ACCH2', 'ACCH', 'C#C', 'CH=CH',
    'CH2=C', 'CH=C', 'C=C', 'OH', 'SiH3', 'CH3', 'SiH2', 'SiH', 'Si',
    'ACH', 'AC', 'CH2', 'CH', 'C'
    ]

# get the fragmentation scheme in the format necessary
fragmentation_scheme = {i+1: j[1] for i, j in enumerate(UNIFAC_SMARTS)}

temp = {j[0]: i+1 for i, j in enumerate(UNIFAC_SMARTS)}
sorted_group_numbers_as_in_paper = [temp[group_name] for group_name in sorted_group_names_as_in_paper]

def function_to_choose_fragmentation(fragmentations):
    fragmentations_descriptors = {}
    i = 0
    for fragmentation in fragmentations:
        fragmentations_descriptors[i] = [len(fragmentation)]
        i += 1
    
    sorted_fragmentations_dict = sorted(fragmentations_descriptors.items(), key=operator.itemgetter(1))

    return fragmentations[sorted_fragmentations_dict[0][0]]

smiles = ['CCCCO', 'CCCO', 'CCO', 'CO']
fragmentation_scheme = {
    'CH2' : '[CH2]',
    'OH' : '[OH]',
    'CH3' : '[CH3]',
    'CH2-CH2' : '[CH2][CH2]'
}
fragmentation_scheme_order1 = ['CH2-CH2', 'CH3', 'CH2', 'OH']

print('simple algorithm 1')
frg = fragmenter(fragmentation_scheme, fragmentation_scheme_order=fragmentation_scheme_order1, algorithm='simple')
for smi in smiles:
    fragmentation, success, fragmentation_matches = frg.fragment(smi)
    print(smi, fragmentation)

# examples of fragmentation drawing
for i, SMILES in enumerate(smiles):
    mol = Chem.MolFromSmiles(SMILES)
    fragmentation, success, fragmentation_matches = frg.fragment(mol)
    img = draw_mol_with_highlights_and_legend(mol, fragmentation_matches)
    img.save(f'simple_example{i+1}.png')

print()
print('simple algorithm 2')
fragmentation_scheme_order2 = ['CH3', 'CH2', 'CH2-CH2', 'OH']
frg = fragmenter(fragmentation_scheme, fragmentation_scheme_order=fragmentation_scheme_order2, algorithm='simple')
for smi in smiles:
    fragmentation, success, fragmentation_matches = frg.fragment(smi)
    print(smi, fragmentation)

print()
print('complete algorithm 1')
frg = fragmenter(fragmentation_scheme, algorithm='complete', n_atoms_cuttoff=30, function_to_choose_fragmentation=function_to_choose_fragmentation)
for smi in smiles:
    fragmentation, success, fragmentation_matches = frg.fragment(smi)
    print(smi, fragmentation)

print()
print('complete algorithm 2')
frg = fragmenter(fragmentation_scheme, algorithm='complete', n_atoms_cuttoff=30, function_to_choose_fragmentation=lambda x: x)
for smi in smiles:
    fragmentations, success, fragmentations_matches = frg.fragment_complete(smi)
    print(smi)
    pprint(fragmentations)
    pprint(fragmentations_matches) # some of the fragmentations are the same, but the found fragmentation_matches are different.