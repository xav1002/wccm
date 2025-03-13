import requests
import re
from io import StringIO
import time
import equilibrator_api as eq
import cobra
import math

def KPA_parse_KEGG(query_items:list[str],cc:eq.ComponentContribution) -> any:
    query_items_all = []
    query_type = query_items[0][0]
    print(query_type)
    for x in query_items:
        if x[0] == 'C' or x[0] == 'R':
            query_items_all.append(x)
        else:
            query_type = 'E'
            query_items_all.append(x)

    metabolites: list = []
    reactions: list = []
    enzymes: list = []

    # STARTHERE: need to break up queries into maximum of 10 at a time
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

        # print(req_raw.text)
        req_2 = list(filter(lambda x: x!='' and x!='\n',re.split(r'///',req_raw.text)))

        metabolites: list = []
        reactions: list = []
        enzymes: list = []

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
                        elif line_level_1.startswith("EXACT_MASS"):
                            category = "EXACT_MASS"
                        elif line_level_1.startswith("REACTION"):
                            category = "REACTION"
                        elif line_level_1.startswith("ENZYME"):
                            category = ""
                        elif line_level_1.startswith("PATHWAY"):
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
                                text = line_level_1.replace("COMMENT","").strip()
                                # In WCCM, need to create map for each generic compound to a list of compounds that fall under that category
                                # Here, just assume that generic acyl-CoA is a hexanoyl-CoA
                                if text == "Generic compound in reaction hierarchy" and entry == 'C00040':
                                    entry = 'C05270'
                                    names = ['Hexanoyl-CoA']
                                    formula = 'C27H46N7O17P3S'
                                    MW = float(865.68)
                                    rxn_names = ['R04747','R04751','R06985','R08796','R08797','R10171']
                                elif text == "Generic compound in reaction hierarchy" and entry == 'C02403':
                                    entry = 'C01585'
                                    names = ['Hexanoic acid','Hexanoate','Hexylic acid','n-Caproic acid']
                                    formula = 'C6H12O2'
                                    MW = float(116.16)
                                    rxn_names = ['R03620']
                    metabolites.append(KPA_Metabolite(entry,names,formula,MW,rxn_names))

                # KPA_Reaction
                case 'R':
                    names: list = []
                    enzyme_id = []

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
                        elif line_level_1.startswith("RCLASS"):
                            category = "RCLASS"
                        elif line_level_1.startswith("ENZYME"):
                            category = "ENZYME"
                        elif line_level_1.startswith("PATHWAY"):
                            category = ""

                        if not line_level_1.startswith("/"):
                            if category == "ENTRY":
                                entry = line_level_1.replace("ENTRY","").replace("Reaction","").strip()
                                # print(entry)
                            elif category == "NAME":
                                names = names+list(filter(lambda x: x!='',list(map(lambda x: x.strip(),re.split(';',line_level_1.replace("NAME",""))))))
                            elif category == "DEFINITION":
                                definition = line_level_1.replace("DEFINITION","").strip()
                            elif line_level_1.startswith("EQUATION"):
                                equation = line_level_1.replace("EQUATION","").strip().replace("C00040","C05270").replace("C02403","C01585")
                            elif line_level_1.startswith("ENZYME"):
                                enzyme_id = list(set(re.split(' ',line_level_1.replace("ENZYME","").strip())))
                                while '' in enzyme_id: enzyme_id.remove('')
                                if '2.3.1.16' in enzyme_id: enzyme_id.remove('2.3.1.16')
                    reactions.append(KPA_Reaction(entry,names,definition,equation,enzyme_id,cc))

                # KPA_Enzyme
                case 'E':
                    names: list = []
                    substrates: list = []
                    products: list = []

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

                        if not line_level_1.startswith("/"):
                            if category == "ENTRY":
                                entry = line_level_1.replace("ENTRY","").replace("Enzyme","").replace("EC ","").strip()
                            elif category == "NAME":
                                names = list(map(lambda x: x.strip(),re.split(';',line_level_1.replace("NAME",""))))
                            elif category == "SYSNAME":
                                sysname = line_level_1.replace("SYSNAME","").strip()
                            elif category == "REACTION":
                                print('test2',line_level_1)
                                reaction_entries = line_level_1.replace(';','').replace("REACTION","").replace("REACTION(KEGG)","").strip()
                                print('test3',line_level_1)
                                reaction_entries = re.split(' ',reaction_entries)
                                while '' in reaction_entries: reaction_entries.remove('')
                                print('test',reaction_entries)
                                for rxn in reaction_entries:
                                    if rxn[0] != 'R':
                                        reaction_entries.remove(rxn)
                            elif category == "SUBSTRATE":
                                substrates = substrates+list(filter(lambda x: x!='',re.split(';',line_level_1.replace("SUBSTRATE","").strip())))
                            elif category == "PRODUCT":
                                products = products+list(filter(lambda x: x!='',re.split(';',line_level_1.replace("PRODUCT","").strip())))
                    enzymes.append(KPA_Enzyme(entry,names,sysname,reaction_entries,substrates,products))

    # attrs = vars(enzymes[0])
    # print(', '.join("%s: %s" % item for item in attrs.items()))

    return (metabolites,reactions,enzymes)

class KPA_Metabolite:

    def __init__(self,entry,names,formula,MW,rxn_names):
        self.__entry = entry
        self.__names = names
        self.__formula = formula
        self.__MW = MW
        self.__rxns = rxn_names

    @property
    def entry(self) -> str:
        return self.__entry

    @property
    def names(self) -> list[str]:
        return self.__names

    @property
    def formula(self) -> str:
        return self.__formula

    @property
    def MW(self) -> float:
        return self.__MW

    @property
    def rxns(self) -> list[str]:
        return self.__rxns

class KPA_Reaction:

    def __init__(self,entry,names,definition,equation,enzyme_id,cc:eq.ComponentContribution):
        self.__entry = entry
        self.__names = names
        self.__definition = definition
        self.__equation = equation
        self.__enzyme_id = enzyme_id

        [self.__stoich,self.__equil_rxns] = self.parse_equation(cc)
        self.__metabolites = []
        self.__cobra_rxn = cobra.Reaction(self.__entry,lower_bound=-100,upper_bound=100,name=self.__enzyme_id[0])
        self.__cobra_rxn.add_metabolites(self.__stoich)

    @property
    def entry(self) -> str:
        return self.__entry

    @property
    def names(self) -> list[str]:
        return self.__names

    @property
    def definition(self) -> str:
        return self.__definition

    @property
    def equation(self) -> str:
        return self.__equation

    @property
    def enzyme_id(self) -> str:
        return self.__enzyme_id

    @property
    def stoich(self) -> dict[str,float]:
        return self.__stoich

    @property
    def equil_rxns(self) -> dict:
        return self.__equil_rxns

    @property
    def metabolites(self) -> list[KPA_Metabolite]:
        return self.__metabolites

    @metabolites.setter
    def metabolites(self,new_metabolites:list[KPA_Metabolite]) -> None:
        self.__metabolites = new_metabolites

    @property
    def cobra_rxn(self) -> cobra.Reaction:
        return self.__cobra_rxn

    def parse_equation(self,cc:eq.ComponentContribution) -> dict:
        [sub_str,prod_str] = re.split('<=>',self.__equation)
        [sub_names,prod_names] = re.split('<=>',self.__definition)

        subs_wth_stoich = re.split(' \\+ ',sub_str)
        prod_wth_stoich = re.split(' \\+ ',prod_str)
        sub_names_wth_stoich = re.split(' \\+ ',sub_names)
        prod_names_wth_stoich = re.split(' \\+ ',prod_names)

        new_stoich = {}
        new_rxn = {}

        for idx,sub in enumerate(subs_wth_stoich):
            stoich = re.split(' ',sub.strip())
            sub_name_stoich = re.split(' ',sub_names_wth_stoich[idx].strip())
            if len(stoich) == 1:
                metabolite = cobra.Metabolite(id=str(stoich[0]),name=sub_name_stoich[0],
                                              compartment='e')
                compound = cc.get_compound(compound_id='kegg:'+str(stoich[0]))
                num = 1
            else:
                name = sub_name_stoich[1]
                if len(sub_name_stoich) == 2:
                    name = sub_name_stoich[1]
                else:
                    name = ''.join(str(x) for x in sub_name_stoich[1:])
                metabolite = cobra.Metabolite(id=str(stoich[1]),name=name,
                                            compartment='e')
                compound = cc.get_compound(compound_id='kegg:'+str(stoich[1]))
                num = stoich[0].replace('n','')
                if num == '': num = '1'
            new_stoich[metabolite] = -int(num)
            new_rxn[compound] = -int(num)
        for idx,prod in enumerate(prod_wth_stoich):
            stoich = re.split(' ',prod.strip())
            prod_name_stoich = re.split(' ',prod_names_wth_stoich[idx].strip())
            if len(stoich) == 1:
                metabolite = cobra.Metabolite(id=str(stoich[0]),name=prod_name_stoich[0],
                                              compartment='e')
                compound = cc.get_compound(compound_id='kegg:'+str(stoich[0]))
                num = 1
            else:
                name = prod_name_stoich[1]
                if len(prod_name_stoich) == 2:
                    name = prod_name_stoich[1]
                else:
                    name = ''.join(str(x) for x in prod_name_stoich[1:])
                metabolite = cobra.Metabolite(id=str(stoich[1]),name=name,
                                              compartment='e')
                compound = cc.get_compound(compound_id='kegg:'+str(stoich[1]))
                num = stoich[0].replace('n','')
                if num == '': num = '1'
            new_stoich[metabolite] = int(num)
            new_rxn[compound] = int(num)
            # print('test6',prod_name_stoich,stoich)

        return [new_stoich,eq.Reaction(new_rxn)]

    def get_metabolite_entries(self) -> list[str]:
        names: list[str] = []
        for meta in list(self.__stoich.keys()):
            names.append(meta.id)
        return names

class KPA_Enzyme:

    def __init__(self,entry,names,sysname,reactions,substrates,products):
        self.__entry = entry
        self.__names = names
        self.__sysname = sysname
        self.__reactions = reactions
        self.__substrates = substrates
        self.__products = products

    @property
    def entry(self) -> str:
        return self.__entry

    @entry.setter
    def entry(self,new_entry) -> None:
        self.__entry = new_entry

    @property
    def names(self) -> list[str]:
        return self.__names

    @property
    def sysname(self) -> str:
        return self.__sysname

    @property
    def reactions(self) -> list[str]:
        return self.__reactions

    @property
    def substrates(self) -> list[str]:
        return self.__substrates

    @property
    def products(self) -> list[str]:
        return self.__products
