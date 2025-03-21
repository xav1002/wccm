import re

def parse_inchi(inchi:str) -> dict[str,int]:
    parsed_inchi = {
        'C_ct':0,
        'H_ct':0,
        'N_ct':0,
        'O_ct':0,
        'P_ct':0,
        'tot_C_ox':0
    }
    layers = {
        'formula':True,
        'c':True,
        'h':True
    }
    formula_str = ''
    c_str = ''
    h_str = ''

    # parse layers
    for key in list(layers.keys()):
        match key:
            case 'formula':
                split_str = inchi.split('/')
                formula_str = split_str[1]
            case 'c':
                split_str = inchi.split('/c')[1].split('/')
                c_str = split_str[0]
            case 'h':
                split_str = inchi.split('/h')[1].split('/')
                h_str = split_str[0]

    print(formula_str,c_str,h_str)

    # count atoms
    letters = re.findall(r'[A-Z]',formula_str)
    nums = re.findall(r'\d+', formula_str)
    if len(letters) > len(nums):
        for idx,x in enumerate(letters):
            letter_idx = formula_str.find(x)
            if letter_idx+1 < len(formula_str) and not formula_str[letter_idx+1].isnumeric():
                nums.insert(idx,'1')
        if not formula_str[-1].isnumeric():
            nums.append('1')
    print(letters,nums)
    for idx,x in enumerate(letters):
        parsed_inchi[x+'_ct'] = int(nums[idx])

    # tally redox of each carbon atom
    for idx,x in enumerate(parsed_inchi['C_ct']):
        
        parsed_inchi['tot_C_ox'] += 1

    return parsed_inchi

# testing parse_inchi
inchi_strs = ['InChI=1S/C4H6O5/c5-2(4(8)9)1-3(6)7/h2,5H,1H2,(H,6,7)(H,8,9)','InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3','InChI=1S/C6H8O6/c7-1-2(8)5-3(9)4(10)6(11)12-5/h2,5,7-10H,1H2/t2-,5+/m0/s1']
smiles_strs = ['C(C(C(=O)O)O)C(=O)O']

import re

def smi_tokenizer(smi):
    """
    Tokenize a SMILES molecule or reaction
    """
    pattern =  "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|=|#|-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"
    regex = re.compile(pattern)
    tokens = [token for token in regex.findall(smi)]
    assert smi == ''.join(tokens)
    return tokens

for x in smiles_strs:
    print(smi_tokenizer(x))