import operator
from fragmenter import fragmenter
from rdkit import Chem
from tqdm import tqdm
import SMARTS

UNIFAC_SMARTS = SMARTS.UNIFAC.copy()

def info_to_CSV(inchikey, SMILES, pubchem_id, fragmentation):
    
    fragmentation_array = []
    for group_number, amount in fragmentation.items():
        fragmentation_array.append(str(group_number) + ":" + str(amount))
    
    return inchikey + "," + SMILES + "," + pubchem_id + "," + "|".join(fragmentation_array)

def CSV_to_info(CSV_line, has_fragmentation = False):
    CSV_line = CSV_line.replace('\n', '')
    array = CSV_line.split(',')
    
    fragmentation = {}
    
    if has_fragmentation:
        fragmentation_array = array[3].split('|')
        for match_str in fragmentation_array:
            array2 = match_str.split(':')
            group_number = int(array2[0])
            amount = int(array2[1])
            
            fragmentation[group_number] = amount
        
    return array[0], array[1], array[2], fragmentation
    
def function_to_choose_fragmentation(fragmentations):
    fragmentations_descriptors = {}
    i = 0
    for fragmentation in fragmentations:
        fragmentations_descriptors[i] = [len(fragmentation)]
        i += 1
    
    sorted_fragmentations_dict = sorted(fragmentations_descriptors.items(), key=operator.itemgetter(1))

    return fragmentations[sorted_fragmentations_dict[0][0]]

def is_fragmentation_equal_to_other_fragmentation(fragmentation, other_fragmentation):
    
    for group_number, amount in fragmentation.items():
        if group_number in other_fragmentation:
            if fragmentation[group_number] != other_fragmentation[group_number]:
                return False
    return True

def log_structure_results(f, pubchem_id, SMILES, inchikey, success, fragmentation, fragmentation_reference_DB, status = ''):
    
    f.write('https://pubchem.ncbi.nlm.nih.gov/compound/' + pubchem_id + '#section=2D-Structure\n')
    f.write(SMILES + '\n')
    f.write(inchikey + '\n')
    f.write('\n' + 'Fragmentation was successfull: ' + str(success) + '\n')
    
    if status != '':
        f.write(status + '\n')
    
    if success:
        f.write('Fragmentation from the algorithm:\n')
        sorted_group_number = sorted(fragmentation.keys())
        
        for group_number in sorted_group_number:
            f.write((UNIFAC_SMARTS[group_number - 1][0]).ljust(12, ' ') + '\t' + str(group_number).ljust(8, ' ') + str(fragmentation[group_number]).ljust(8, ' ') + '\n')
    
    f.write('\n')
    
    if len(fragmentation_reference_DB) > 0:
        f.write('Fragmentation from the reference database:\n')
        sorted_group_number = sorted(fragmentation_reference_DB.keys())
        
        for group_number in sorted_group_number:
            f.write((UNIFAC_SMARTS[group_number - 1][0]).ljust(12, ' ') + '\t' + str(group_number).ljust(8, ' ') + str(fragmentation_reference_DB[group_number]).ljust(8, ' ') + '\n')
    
        
    f.write('\n\n')   


# the groups in groups_not_on_DDB_list are not available
# on the list of the DDB: https://www.ddbst.com/published-parameters-unifac.html
# This is why they are taken out of the fragmentation scheme, to reproduce the results
# of the reference database. In case you don't care about reproducing the results from the paper
# but being able to fragent more molecules, you can empty the list, so that these groups are 
# included in the fragmentation scheme or have a look at execute_fragmentation.py
groups_not_on_DDB_list = ["H2COCH", "HCOCH", "COCH", "H2COCH2", "OCOCO", "(CH3O)2CO", "(CH2O)2CO", "CH3OCH2OCO"]

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
for group_name in groups_not_on_DDB_list:
    sorted_group_names_as_in_paper.remove(group_name)
    for i, (gn, sm) in enumerate(UNIFAC_SMARTS):
        if gn == group_name:
            UNIFAC_SMARTS[i] = (group_name, "")

# get the fragmentation scheme in the format necessary
fragmentation_scheme = {i+1: j[1] for i, j in enumerate(UNIFAC_SMARTS)}

temp = {j[0]: i+1 for i, j in enumerate(UNIFAC_SMARTS)}
sorted_group_numbers_as_in_paper = [temp[group_name] for group_name in sorted_group_names_as_in_paper]

if set(sorted_group_names_as_in_paper).union(set(groups_not_on_DDB_list)) != set([t[0] for t in UNIFAC_SMARTS]) \
        or len(sorted_group_names_as_in_paper) + len(groups_not_on_DDB_list) != len(UNIFAC_SMARTS):
    raise ValueError('Missmatch of groups. If you are getting this error, \
                     you probably changed a group name, deleted or added one group. \
                     To avoid this error reflect the changes you made in "sorted_group_names_as_in_paper"')

# first step: fragment reference database and compare with the results
reference_DB = []
with open('reference_DB.csv') as f:
    for line in f.readlines():
        reference_DB.append(CSV_to_info(line, True))

reference_DB_fragmentation_stats = {}

simple_fragmenter_fragmented = []
simple_fragmenter_fragmented_and_equal_to_reference_DB = []
complete_fragmenter_fragmented = []
complete_fragmenter_fragmented_and_equal_to_reference_DB = []

right_size_for_complete_fragmenter = []

simple_fragmenter = fragmenter(fragmentation_scheme, algorithm='simple')
complete_fragmenter = fragmenter(fragmentation_scheme, algorithm='complete', n_atoms_cuttoff=20, function_to_choose_fragmentation=function_to_choose_fragmentation)

# without sorting the patterns
print('####################################################################')
print('Fragmenting the reference database without the patterns sorted (simple and complete algorithm)')

# I don't know why, but for the unordered group search I am not getting the same result as the old algorithm
with (open('reference_DB_simple_fragmentation_without_pattern_sorting_results.log','w+') as f_simple,
      open('reference_DB_complete_fragmentation_without_pattern_sorting_results.log','w+') as f_complete):
    for inchikey, SMILES, pubchem_id, fragmentation_reference_DB in tqdm(reference_DB, total=len(reference_DB)):
        lines = []
        
        
        for group_number, amount in fragmentation_reference_DB.items():
            if not group_number in reference_DB_fragmentation_stats:
                reference_DB_fragmentation_stats[group_number] = 0
                
            reference_DB_fragmentation_stats[group_number] += amount
        
        fragmentation, success, fragmentation_matches = simple_fragmenter.fragment(SMILES)
        if success:
            simple_fragmenter_fragmented.append(inchikey)
            if is_fragmentation_equal_to_other_fragmentation(fragmentation, fragmentation_reference_DB):
                simple_fragmenter_fragmented_and_equal_to_reference_DB.append(inchikey)
        
        log_structure_results(f_simple, pubchem_id, SMILES, inchikey, success, fragmentation, fragmentation_reference_DB)
        
        n_heavy_atoms  = 0
        for sub_SMILES in SMILES.split("."):
            n_heavy_atoms = max(n_heavy_atoms, fragmenter.get_heavy_atom_count(Chem.MolFromSmiles(sub_SMILES)))
        
        if n_heavy_atoms <= 20:
            right_size_for_complete_fragmenter.append(inchikey)
            fragmentation, success, fragmentation_matches = complete_fragmenter.fragment(SMILES)
            if success:
                complete_fragmenter_fragmented.append(inchikey)
                if is_fragmentation_equal_to_other_fragmentation(fragmentation, fragmentation_reference_DB):
                    complete_fragmenter_fragmented_and_equal_to_reference_DB.append(inchikey)
            
            log_structure_results(f_complete, pubchem_id, SMILES, inchikey, success, fragmentation, fragmentation_reference_DB)
        else:
            log_structure_results(f_complete, pubchem_id, SMILES, inchikey, success, fragmentation, fragmentation_reference_DB, 'Structure was skipped because it is larger than 20 atoms.')


print('')
print('N_structures(simple): ' + str(len(reference_DB)))
print('N_fragmented(simple): ' + str(len(simple_fragmenter_fragmented)) + "(" + str((1.0 * len(simple_fragmenter_fragmented)) / len(reference_DB)) + ")")
print('N_fragmented_and_equal(simple): ' + str(len(simple_fragmenter_fragmented_and_equal_to_reference_DB)) + "(" + str((1.0 * len(simple_fragmenter_fragmented_and_equal_to_reference_DB)) / len(reference_DB)) + ")")
print('')
print('N_structures(complete):' + str(len(right_size_for_complete_fragmenter)))
print('N_fragmented(complete): ' + str(len(complete_fragmenter_fragmented)) + "(" + str((1.0 * len(complete_fragmenter_fragmented)) / len(right_size_for_complete_fragmenter)) + ")")
print('N_fragmented_and_equal(complete): ' + str(len(complete_fragmenter_fragmented_and_equal_to_reference_DB)) + "(" + str((1.0 * len(complete_fragmenter_fragmented_and_equal_to_reference_DB)) / len(right_size_for_complete_fragmenter)) + ")")
print('')
print('')
print('')
print('')

right_size_for_complete_fragmenter2 = []

# with sorting the patterns
simple_fragmenter.fragmentation_scheme_order = sorted_group_numbers_as_in_paper
complete_fragmenter.fragmentation_scheme_order = sorted_group_numbers_as_in_paper

simple_fragmenter_sorted_fragmented = []
simple_fragmenter_sorted_fragmented_and_equal_to_reference_DB = []
complete_fragmenter_sorted_fragmented = []
complete_fragmenter_sorted_fragmented_and_equal_to_reference_DB = []


print('####################################################################')
print('Fragmenting the reference database with the patterns sorted (simple and complete algorithm)')
# with the sorted groups I am getting the same result as the old algorithm
with (open('reference_DB_simple_fragmentation_with_pattern_sorting_results.log','w+') as f_simple,
      open('reference_DB_complete_fragmentation_with_pattern_sorting_results.log','w+') as f_complete):
    for inchikey, SMILES, pubchem_id, fragmentation_reference_DB in tqdm(reference_DB, total=len(reference_DB)):
        
        fragmentation, success, fragmentation_matches = simple_fragmenter.fragment(SMILES)
        if success:
            simple_fragmenter_sorted_fragmented.append(inchikey)
            if is_fragmentation_equal_to_other_fragmentation(fragmentation, fragmentation_reference_DB):
                simple_fragmenter_sorted_fragmented_and_equal_to_reference_DB.append(inchikey)
                
        log_structure_results(f_simple, pubchem_id, SMILES, inchikey, success, fragmentation, fragmentation_reference_DB)
        
        n_heavy_atoms  = 0
        for sub_SMILES in SMILES.split("."):
            n_heavy_atoms = max(n_heavy_atoms, fragmenter.get_heavy_atom_count(Chem.MolFromSmiles(sub_SMILES)))
        
        if n_heavy_atoms <= 20:
            right_size_for_complete_fragmenter2.append(inchikey)
            fragmentation, success, fragmentation_matches = complete_fragmenter.fragment(SMILES)
            if success:
                complete_fragmenter_sorted_fragmented.append(inchikey)
                if is_fragmentation_equal_to_other_fragmentation(fragmentation, fragmentation_reference_DB):
                    complete_fragmenter_sorted_fragmented_and_equal_to_reference_DB.append(inchikey)
            
            log_structure_results(f_complete, pubchem_id, SMILES, inchikey, success, fragmentation, fragmentation_reference_DB)
        else:
            log_structure_results(f_complete, pubchem_id, SMILES, inchikey, success, fragmentation, fragmentation_reference_DB, 'Structure was skipped because it is larger than 20 atoms.')
 
print('')
print('N_structures(simple): ' + str(len(reference_DB)))
print('N_fragmented(simple): ' + str(len(simple_fragmenter_sorted_fragmented)) + "(" + str((1.0 * len(simple_fragmenter_sorted_fragmented)) / len(reference_DB)) + ")")
print('N_fragmented_and_equal(simple): ' + str(len(simple_fragmenter_sorted_fragmented_and_equal_to_reference_DB)) + "(" + str((1.0 * len(simple_fragmenter_sorted_fragmented_and_equal_to_reference_DB)) / len(reference_DB)) + ")")
print('')
print('N_structures(complete):' + str(len(right_size_for_complete_fragmenter2)))
print('N_fragmented(complete): ' + str(len(complete_fragmenter_sorted_fragmented)) + "(" + str((1.0 * len(complete_fragmenter_sorted_fragmented)) / len(right_size_for_complete_fragmenter2)) + ")")
print('N_fragmented_and_equal(complete): ' + str(len(complete_fragmenter_sorted_fragmented_and_equal_to_reference_DB)) + "(" + str((1.0 * len(complete_fragmenter_sorted_fragmented_and_equal_to_reference_DB)) / len(right_size_for_complete_fragmenter2)) + ")")
print('')
print('')
print('')
print('')


# second step: try to fragent all from the component    
structures_DB = []
with open('structures_DB.csv') as f:
    for line in f.readlines():
        structures_DB.append(CSV_to_info(line))
        
combined_fragmenter = fragmenter(fragmentation_scheme, algorithm='combined', n_atoms_cuttoff=20, function_to_choose_fragmentation=function_to_choose_fragmentation, n_max_fragmentations_to_find=1)

combined_fragmenter.fragmentation_scheme_order = sorted_group_numbers_as_in_paper 
combined_fragmenter.n_max_fragmentations_to_find = 1    

combined_fragmenter_sorted_fragmented = []
right_size_for_combined_fragmenter  = []
print('####################################################################')
print('Fragmenting the structures database with the patterns sorted (combined algorithm)')
with open('structures_DB_combined_fragmentation_with_pattern_sorting_results.log','w+') as f_combined:
    for inchikey, SMILES, pubchem_id, empty_fragmentation in tqdm(structures_DB, total=len(structures_DB)):
        
        fragmentation, success, fragmentation_matches = combined_fragmenter.fragment(SMILES)
        if success:
            combined_fragmenter_sorted_fragmented.append(inchikey)
            
        n_heavy_atoms  = 0
        for sub_SMILES in SMILES.split("."):
            n_heavy_atoms = max(n_heavy_atoms, fragmenter.get_heavy_atom_count(Chem.MolFromSmiles(sub_SMILES)))
        
        if n_heavy_atoms <= 20:
            right_size_for_combined_fragmenter.append(inchikey)
            log_structure_results(f_combined, pubchem_id, SMILES, inchikey, success, fragmentation, {})
        else:
            if success:
                log_structure_results(f_combined, pubchem_id, SMILES, inchikey, success, fragmentation, {})
            else:
                log_structure_results(f_combined, pubchem_id, SMILES, inchikey, success, fragmentation, {}, 'Structure was skipped because it is larger than 20 atoms.')
    
print('')
print('N_structures(simple): ' + str(len(structures_DB)))
print('N_fragmented(simple): ' + str(len(combined_fragmenter_sorted_fragmented)) + "(" + str((1.0 * len(combined_fragmenter_sorted_fragmented)) / len(structures_DB)) + ")")
print('')
print('####################################################################')