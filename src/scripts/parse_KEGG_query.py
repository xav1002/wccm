import requests
import re
from io import StringIO
import time
from equilibrator_api import ComponentContribution
import math

from src.classes.components.MPNG_Metabolite import MPNG_Metabolite
from src.classes.components.MPNG_Reaction import MPNG_Reaction
from src.classes.components.MPNG_Enzyme import MPNG_Enzyme

def parse_KEGG(query_items:list[str]) -> MPNG_Metabolite | MPNG_Reaction | MPNG_Enzyme:
    query_items_all = []
    query_type = query_items[0][0]
    for x in query_items:
        if x[0] == 'C' or x[0] == 'R':
            query_items_all.append(x)
        else:
            query_type = 'E'
            query_items_all.append(x)

    metabolites: list = []
    reactions: list = []
    enzymes: list = []

    for idx in range(math.ceil(len(query_items_all)/10)):
        if len(query_items_all) < 10:
            query_items_segment = query_items_all
        else:
            query_items_segment = query_items_all[0:10]
            query_items_all = [item for item in query_items_all if item not in query_items_segment]

        query_items_str = query_items_segment.pop(0)
        for item in query_items_segment:
            query_items_str = query_items_str+'+'+item

        KEGG_link = 'https://rest.kegg.jp/get/'+query_items_str

        print(KEGG_link)
        try:
            req_raw = requests.get(KEGG_link)
        except Exception as e:
            print('KEGG request error:',e)
            time.sleep(1/3)
            req_raw = requests.get(KEGG_link)

        req_2 = list(filter(lambda x: x!='' and x!='\n',re.split(r'///',req_raw.text)))

        for x in req_2:
            req_3 = StringIO(x)
            category = ""

            match query_type:
                # Compound
                case 'C':
                    names: list = []
                    rxn_names: list = []
                    MW = 0.0
                    formula = ''
                    is_generic = False
                    BRITE_lvl = 0
                    BRITE_dict = {}
                    DBLINKS_dict = {}
                    for line_level_1 in req_3:
                        line_level_1 = line_level_1.strip()

                        if line_level_1.startswith("ENTRY"):
                            category = "ENTRY"
                        elif line_level_1.startswith("NAME"):
                            category = "NAME"
                        elif line_level_1.startswith("FORMULA"):
                            category = "FORMULA"
                        elif line_level_1.startswith("MOL_WEIGHT"):
                            category = "MOL_WEIGHT"
                        elif line_level_1.startswith("REMARK"):
                            category = ""
                        elif line_level_1.startswith("COMMENT"):
                            category = "COMMENT"
                        elif line_level_1.startswith("MODULE"):
                            category = ""
                        elif line_level_1.startswith("EXACT_MASS"):
                            category = "EXACT_MASS"
                        elif line_level_1.startswith("REACTION"):
                            category = "REACTION"
                        elif line_level_1.startswith("ENZYME"):
                            category = ""
                        elif line_level_1.startswith("PATHWAY"):
                            category = ""
                        elif line_level_1.startswith("SEQUENCE"):
                            category = ""
                        elif line_level_1.startswith("GENE"):
                            category = ""
                        elif line_level_1.startswith("ORGANISM"):
                            category = ""
                        elif line_level_1.startswith("TYPE"):
                            category = ""
                        elif line_level_1.startswith("BRITE"):
                            category = "BRITE"
                        elif line_level_1.startswith("ATOM"):
                            category = ""
                        elif line_level_1.startswith("BOND"):
                            category = ""
                        elif line_level_1.startswith("PATHWAY"):
                            category = ""
                        elif line_level_1.startswith("DBLINKS"):
                            category = "DBLINKS"
                        elif line_level_1.startswith("NETWORK"):
                            category = ""

                        if not line_level_1.startswith("/"):
                            if category == "ENTRY":
                                entry = line_level_1.replace("ENTRY","").replace("Compound","").strip()
                            elif category == "NAME":
                                names = names+list(filter(lambda x: x!='',list(map(lambda x: x.strip(),re.split(';',line_level_1.replace("NAME",""))))))
                            elif category == "FORMULA":
                                formula = line_level_1.replace("FORMULA","").strip()
                            elif category == "MOL_WEIGHT":
                                MW = float(line_level_1.replace("MOL_WEIGHT","").strip())
                            elif category == "REACTION":
                                rxn_names = rxn_names+list(filter(lambda x: x!='',list(map(lambda x: x.strip(),re.split(' ',line_level_1.replace("REACTION",""))))))
                            elif category == "COMMENT":
                                if 'Generic compound in reaction hierarchy' in line_level_1:
                                    is_generic = True
                                else:
                                    is_generic = False
                            elif category == "BRITE":
                                BRITE_lvl += 1
                                BRITE_dict[BRITE_lvl] = line_level_1.replace("BRITE","").strip()
                            elif category == "DBLINKS":
                                line_level_1 = line_level_1.replace("DBLINKS","").strip()
                                DBLINKS_dict[line_level_1.split(":")[0]] = line_level_1.split(":")[1].strip()
                    metabolites.append(MPNG_Metabolite(entry,names,formula,MW,rxn_names,is_generic,BRITE_dict))

                # MPNG_Reaction
                case 'R':
                    names: list = []
                    enzyme_id = ['']

                    for line_level_1 in req_3:
                        line_level_1 = line_level_1.strip()

                        if line_level_1.startswith("ENTRY"):
                            category = "ENTRY"
                        elif line_level_1.startswith("NAME"):
                            category = "NAME"
                        elif line_level_1.startswith("DEFINITION"):
                            category = "DEFINITION"
                        elif line_level_1.startswith("EQUATION"):
                            category = "EQUATION"
                        elif line_level_1.startswith("COMMENT"):
                            category = ""
                        elif line_level_1.startswith("RCLASS"):
                            category = "RCLASS"
                        elif line_level_1.startswith("ENZYME"):
                            category = "ENZYME"
                        elif line_level_1.startswith("PATHWAY"):
                            category = ""
                        elif line_level_1.startswith("DBLINKS"):
                            category = ""
                        elif line_level_1.startswith("NETWORK"):
                            category = ""

                        if not line_level_1.startswith("/"):
                            if category == "ENTRY":
                                entry = line_level_1.replace("ENTRY","").replace("Reaction","").replace("Overall","").strip()
                            elif category == "NAME":
                                names = names+list(filter(lambda x: x!='',list(map(lambda x: x.strip(),re.split(';',line_level_1.replace("NAME",""))))))
                            elif category == "DEFINITION":
                                definition = line_level_1.replace("DEFINITION","").strip()
                            elif line_level_1.startswith("EQUATION"):
                                equation = line_level_1.replace("EQUATION","").strip()
                            elif line_level_1.startswith("ENZYME"):
                                enzyme_id = list(set(re.split(' ',line_level_1.replace("ENZYME","").strip())))
                                while '' in enzyme_id: enzyme_id.remove('')
                    reactions.append(MPNG_Reaction(entry,names,definition,equation,enzyme_id,True,True))

                # MPNG_Enzyme
                case 'E':
                    names: list = []
                    enz_rxn_ids: list = []
                    substrates: list = []
                    products: list = []
                    sysname = ''

                    for line_level_1 in req_3:
                        line_level_1 = line_level_1.strip()

                        if line_level_1.startswith("ENTRY"):
                            category = "ENTRY"
                        elif line_level_1.startswith("NAME"):
                            category = "NAME"
                        elif line_level_1.startswith("CLASS"):
                            category = "CLASS"
                        elif line_level_1.startswith("SYSNAME"):
                            category = "SYSNAME"
                        elif line_level_1.startswith("REACTION"):
                            category = ""
                        elif line_level_1.startswith("ALL_REAC"):
                            category = "REACTION"
                        elif line_level_1.startswith("SUBSTRATE"):
                            category = "SUBSTRATE"
                        elif line_level_1.startswith("PRODUCT"):
                            category = "PRODUCT"
                        elif line_level_1.startswith("COMMENT"):
                            category = ""
                        elif line_level_1.startswith("PATHWAY"):
                            category = ""
                        elif line_level_1.startswith("DBLINKS"):
                            category = ""
                        elif line_level_1.startswith("NETWORK"):
                            category = ""

                        if not line_level_1.startswith("/"):
                            if category == "ENTRY":
                                entry = line_level_1.replace("ENTRY","").replace("Enzyme","").replace("EC ","").strip()
                            elif category == "NAME":
                                names = list(map(lambda x: x.strip(),re.split(';',line_level_1.replace("NAME",""))))
                            elif category == "SYSNAME":
                                sysname = line_level_1.replace("SYSNAME","").strip()
                            elif category == "SUBSTRATE":
                                substrates = substrates+list(filter(lambda x: x!='',re.split(';',line_level_1.replace("SUBSTRATE","").strip())))
                            elif category == "PRODUCT":
                                products = products+list(filter(lambda x: x!='',re.split(';',line_level_1.replace("PRODUCT","").strip())))
                            elif category == "REACTION":
                                enz_rxn_ids = enz_rxn_ids+list(filter(lambda x: x!='',re.split(' ',line_level_1.
                                                                                           replace("ALL_REAC","").
                                                                                           replace("(other)","").
                                                                                           replace(">"," ").
                                                                                           replace(";","").strip())))
                    enzymes.append(MPNG_Enzyme(entry,names,sysname,enz_rxn_ids,substrates,products))
                case 'Glycan':
                    test = 1

        # attrs = vars(enzymes[0])
        # print(', '.join("%s: %s" % item for item in attrs.items()))

    print([metabolites,reactions,enzymes])
    return [metabolites,reactions,enzymes]