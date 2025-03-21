from pyvis.network import Network
import networkx as nx
import json
import pandas as pd
import re
import copy

from cobra import Metabolite

from zeep import Client
import hashlib
import requests
import time
import os
import sys

import zeep
import zeep.helpers

from src.classes.MetabolicPathwayNetworkGraph import MetabolicPathwayNetworkGraph
from src.classes.MicrobialConsortiumKineticModel import MicrobialConsortiumKineticModel

from src.scripts.parse_KEGG_query import parse_KEGG

from src.classes.components.MPNG_Metabolite import MPNG_Metabolite
from src.classes.components.MPNG_Reaction import MPNG_Reaction
from src.classes.components.MPNG_Enzyme import MPNG_Enzyme

class WholeCellConsortiumModel:
    def __init__(self,valid_cofactors:dict[str,bool]):
        self.__metabolites: dict[str,MPNG_Metabolite] = {}
        self.__reactions: dict[str,MPNG_Reaction] = {}
        self.__enzymes: dict[str,MPNG_Enzyme] = {}

        self.__db_excluded_metas: list[str] = ['h+','hydron']
        self.__MPNG_Metabolite_to_CID = {}
        self.__BRENDA_ligand_name_to_CID = {}

        with open('./src/stores/metabolites.json') as f:
            meta_data = json.load(f)
            for x in meta_data:
                new_meta = MPNG_Metabolite.fromJSON(dict_from_json=json.loads(x))
                self.__metabolites[new_meta.entry] = new_meta
        with open('./src/stores/reactions.json') as f:
            rxn_data = json.load(f)
            for x in rxn_data:
                new_rxn = MPNG_Reaction.fromJSON(dict_from_json=json.loads(x))
                self.__reactions[new_rxn.entry] = new_rxn
        with open('./src/stores/enzymes.json') as f:
            enz_data = json.load(f)
            for x in enz_data:
                new_enz = MPNG_Enzyme.fromJSON(dict_from_json=json.loads(x))
                self.__enzymes[new_enz.entry] = new_enz
        with open('./src/stores/KEGG_CIDs.json') as f:
            self.__MPNG_Metabolite_to_CID = json.load(f)
        print('number of metabolites: ',len(list(self.__metabolites.keys())))
        print('number of reactions: ',len(list(self.__reactions.keys())))

        self.__small_gas_metas_entries = ['C00001','C00007','C00011']
        self.__common_metabolite_entries = ['C00138','C00139','C00080','C00024','C00125','C00126']
        for x in range(14):
            zeros = '0'*(5-len(str(x+1)))
            self.__common_metabolite_entries.append('C'+zeros+str(x+1))
        self.common_metabolites = self.__get_metabolites(entries=self.__common_metabolite_entries)

        self.__graphs: dict[str,MetabolicPathwayNetworkGraph] = {}
        self.__whole_KEGG_graph = nx.DiGraph()

        self.__excluded_metabolite_entries = ['C00002','C00008'] + [x for x in list(valid_cofactors.keys()) if not valid_cofactors[x]]
        for meta in list(self.__metabolites.values()):
            if meta.generic:
                self.__excluded_metabolite_entries.append(meta.entry)

        for x in self.__get_reactions('all'):
            meta_ids = list(map(lambda y: y.id,list(x.stoich.keys())))
            if all([z not in meta_ids and 'G' not in z for z in self.__excluded_metabolite_entries]):
                self.add_reaction(x)

        # this is set based on manual curation of KEGG BRITE heirarchy
        self.__generic_compounds = []
        self.__generic_compound_assignments = {

        }

    @property
    def metabolites(self) -> dict[str,MPNG_Metabolite]:
        return self.__metabolites

    @property
    def reactions(self) -> dict[str,MPNG_Reaction]:
        return self.__reactions

    @property
    def enzymes(self) -> dict[str,MPNG_Enzyme]:
        return self.__enzymes

    @metabolites.setter
    def metabolites(self,metas:dict[str,MPNG_Metabolite]) -> None:
        self.__metabolites = metas

    @reactions.setter
    def reactions(self,rxns:dict[str,MPNG_Reaction]) -> None:
        self.__reactions = rxns

    @enzymes.setter
    def enzymes(self,enz:dict[str,MPNG_Enzyme]) -> None:
        self.__enzymes = enz

    @property
    def generic_compound_assignments(self) -> dict[str,str]:
        return self.__generic_compound_assignments

    def __get_metabolites(self,entries:str|list) -> MPNG_Metabolite | list[MPNG_Metabolite]:
        if type(entries) == str:
            if entries == 'all':
                return list(self.__metabolites.values())
            else:
                return self.__metabolites[entries]
        elif type(entries) == list:
            metas = []
            for x in entries:
                metas.append(self.__metabolites[x])
            return metas
        else:
            return []

    def __get_reactions(self,entries:str|list) -> MPNG_Reaction | list[MPNG_Reaction]:
        if type(entries) == str:
            if entries == 'all':
                return list(self.__reactions.values())
            else:
                return self.__reactions[entries]
        elif type(entries) == list:
            rxns = []
            for x in entries:
                rxns.append(self.__reactions[x])
            return rxns
        else:
            return []

    def __get_enzymes(self,entries:str|list) -> MPNG_Enzyme | list[MPNG_Enzyme]:
        if type(entries) == str:
            if entries == 'all':
                return list(self.__enzymes.values())
            else:
                return self.__enzymes[entries]
        elif type(entries) == list:
            enz = []
            for x in entries:
                enz.append(self.__enzymes[x])
            return enz
        else:
            return []

    def add_reaction(self,new_reaction:MPNG_Reaction) -> None:
        # add to NX Graph
        self.__whole_KEGG_graph.add_node(node_for_adding=new_reaction.entry+'_f')
        self.__whole_KEGG_graph.add_node(node_for_adding=new_reaction.entry+'_r')

        stoich: dict[Metabolite,int] = new_reaction.stoich
        key_entries = list(map(lambda x: x.id,list(stoich.keys())))
        for m in key_entries:
            if m not in self.__whole_KEGG_graph.nodes:
                self.__whole_KEGG_graph.add_node(node_for_adding=m)
            if stoich[list(stoich.keys())[key_entries.index(m)]] < 0:
                self.__whole_KEGG_graph.add_edge(m,new_reaction.entry+'_f',
                    stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                )
                self.__whole_KEGG_graph.add_edge(new_reaction.entry+'_r',m,
                    stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                )
            elif stoich[list(stoich.keys())[key_entries.index(m)]] > 0:
                self.__whole_KEGG_graph.add_edge(new_reaction.entry+'_f',m,
                    stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                )
                self.__whole_KEGG_graph.add_edge(m,new_reaction.entry+'_r',
                    stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                )

    def identify_generic_compounds(self) -> pd.DataFrame:
        generic_metas = pd.DataFrame([],columns=['entry','names'])
        for idx,meta in enumerate(list(self.__metabolites.values())):
            if meta.generic:
                generic_metas.loc[idx] = [meta.entry,meta.names]
        print(generic_metas.to_markdown())
        return self.generic_compound_assignments

    def generate_generic_rxns(self) -> None:
        # adding reactions of generic compounds to metas and creating dict of generic-to-specific metas
        generic_to_specific_metas = {}
        for meta in list(self.__metabolites.values()):
            for lvl in list(meta.BRITE_dict.values()):
                if lvl in list(self.__generic_compound_assignments.keys()):
                    if self.__generic_compound_assignments[lvl] not in list(generic_to_specific_metas.keys()):
                        generic_to_specific_metas[self.__generic_compound_assignments[lvl]] = []
                    generic_to_specific_metas[self.__generic_compound_assignments[lvl]] += meta.entry
                    meta.add_reactions(self.__generic_compound_assignments[lvl])

        ### STARTHERE: how to do this?
        # creating MPNG_Reaction objects for specified generic reactions
        generic_rxns = {}
        for meta in list(self.__metabolites.values()):
            for rxn in meta.reactions:
                generic_rxns[rxn] += meta.entry

        for rxn in list(generic_rxns.keys()):
            return
        return

    def match_KEGG_compounds_to_BRENDA(self) -> None:
        # 1. find CID for each KEGG compound from name in KEGG entry
        # 1.1. if no PubChem SID in KEGG entry, try CAS number, if no CAS, then need to use CheBI
        for meta in self.__get_metabolites(['C00256']):
        # for meta in list(self.__metabolites.values()):
            try:
                raw_req = requests.get('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'+meta.names[0]+'/cids/JSON').text
                parsed_req = json.loads(raw_req)
                CIDs = parsed_req['IdentifierList']['CID']
                self.__MPNG_Metabolite_to_CID[meta.entry] = CIDs
                print('meta_names',meta.names[0],CIDs)

                # matching enantiomers to non-steric specified compound
                if any([x in meta.names[0] for x in ['D-','L-','(R)','(S)']]):
                    print('worked')
                    matches = ['(R)-','(S)-','-(R)','-(S)','D-','L-']
                    new_name = meta.names[0]
                    for x in matches:
                        new_name = new_name.replace(x,'')
                    print('new_name',new_name)
                    raw_req_2 = requests.get('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'+new_name+'/cids/JSON').text
                    parsed_req_2 = json.loads(raw_req_2)
                    base_CIDs = parsed_req_2['IdentifierList']['CID']
                    append_to_base = [x for x in list(self.__MPNG_Metabolite_to_CID.keys()) if all([y in self.__MPNG_Metabolite_to_CID[x] for y in base_CIDs])]
                    print('append',append_to_base)
                    for x in append_to_base:
                        self.__MPNG_Metabolite_to_CID[x] += CIDs
                        self.__MPNG_Metabolite_to_CID[x] = list(set(self.__MPNG_Metabolite_to_CID[x]))
                        self.__MPNG_Metabolite_to_CID[meta.entry] += self.__MPNG_Metabolite_to_CID[x]
                        self.__MPNG_Metabolite_to_CID[meta.entry] = list(set(self.__MPNG_Metabolite_to_CID[meta.entry]))
            except Exception as e:
                print('CID not found',e)
        print('MPNG_Metabolite_to_CID',self.__MPNG_Metabolite_to_CID)

        with open('./src/stores/KEGG_CIDs.json', 'w', encoding='utf-8') as f:
            json.dump(self.__MPNG_Metabolite_to_CID, f, ensure_ascii=False, indent=4)

        # NOTE: approach to accumulating the relationships between KEGG and BRENDA compounds/ligands:
        # go through all enzymes, build dict mapping compounds/ligands to their CIDs, save these dicts for later usage
        # for now, ignore generic ligands in BRENDA, assume reversibility
        # NOTE: approach for incorporating generic reaction compounds:
        # make the reversibility the same as the reversibility of a specific compound with respect to a given enzymatic reaction

    def set_reaction_mass_balance(self) -> None:
        ### STARTHERE: need to make exception for generic reaction compounds
        for rxn in self.__get_reactions('all'):
            try:
                stoich = rxn.stoich
                subs = []
                prods = []
                for meta in list(stoich.keys()):
                    if stoich[meta] < 0:
                        subs.append(meta)
                    else:
                        prods.append(meta)
                MPNG_subs = list(map(lambda x: self.__get_metabolites(x.id).MW*abs(stoich[x]),subs))
                MPNG_prods = list(map(lambda x: self.__get_metabolites(x.id).MW*abs(stoich[x]),prods))
                tot_subs_MW = sum(MPNG_subs)
                tot_prods_MW = sum(MPNG_prods)
                if not any([x.generic for x in self.__get_metabolites([y.id for y in subs+prods])]) and (round(tot_subs_MW) != round(tot_prods_MW)):
                    rxn.forward_valid = False
                    rxn.backward_valid = False
            except Exception as e:
                print('mass balance in rxn err, ',e)
                pass

        print('updating reactions.json')
        with open('./src/stores/reactions.json', 'w', encoding='utf-8') as f:
            json.dump(list(map(lambda x: x.toJSON(),list(self.__reactions.values()))), f, ensure_ascii=False, indent=4)

    def set_reaction_reversibility(self) -> None:
        # generating dict mapping MPNG_Metabolites to CID
        # self.match_KEGG_compounds_to_BRENDA()

        # checking whether each enzyme is reversible
        wsdl = "https://www.brenda-enzymes.org/soap/brenda_zeep.wsdl"
        password = hashlib.sha256("b3br?B$iDjpeJm77".encode("utf-8")).hexdigest()
        client = Client(wsdl)
        # for enz in [x for x in list(map(lambda x: x.entry, list(self.__enzymes.values()))) if x[0] == '6']:
        for enz in ['2.3.1.5']:
            # logic:
            # 1. for each reaction that each enzyme catalyzes, are they reversible?
            #   1.a. Track this in MPNG_Reaction
            # 2. If an MPNG_Reaction has >=1 enzyme catalyzing it that is reversible, then the reaction is reversible
            # NOTE: can later extend this to track the organism that an enzyme is from
            # NOTE: this extension will also consider generic reaction compounds. Certain versions of the enzymes can catalyze 
            # reactions to different levels of generalizability, meaning control over the reaction flux directions can be controlled 
            # by using specific enzymes from different organisms.
            try:
                print('starting enz: ',enz)

                # getting natural substrates
                params_nat_sub = ("v.a.xu@wustl.edu",password,f"ecNumber*{enz}", "naturalSubstrate*", "naturalReactionPartners*", "organism*", "ligandStructureId*")
                try:
                    res_nat_sub = client.service.getNaturalSubstrate(*params_nat_sub)
                except Exception as e:
                    print('BRENDA request error:',e)
                    time.sleep(1/3)
                    res_nat_sub = client.service.getNaturalSubstrate(*params_nat_sub)
                # getting substrates
                params_sub = ("v.a.xu@wustl.edu",password,f"ecNumber*{enz}", "substrate*", "reactionPartners*", "organism*", "ligandStructureId*")
                try:
                    res_sub = client.service.getSubstrate(*params_sub)
                except Exception as e:
                    print('BRENDA request error:',e)
                    time.sleep(1/3)
                    res_sub = client.service.getSubstrate(*params_sub)
                # query naturalSubstrate and naturalReactionPartners
                natRxnPartners = [x['naturalReactionPartners'] for x in res_nat_sub] + [x['reactionPartners'] for x in res_sub]
                naturalSubstrate = [x['naturalSubstrate'] for x in res_nat_sub] + [x['substrate'] for x in res_sub]

                # getting natural substrates products
                params_nat_rev = ("v.a.xu@wustl.edu",password,f"ecNumber*{enz}", "organism*", "naturalSubstrates*", "organismNaturalSubstrates*", "commentaryNaturalSubstrates*", "naturalProducts*", "commentaryNaturalProducts*", "organismNaturalProducts*", "reversibility*")
                try:
                    res_nat_rev = client.service.getNaturalSubstratesProducts(*params_nat_rev)
                except Exception as e:
                    print('BRENDA request error:',e)
                    time.sleep(1/3)
                    res_nat_rev = client.service.getNaturalSubstratesProducts(*params_nat_rev)
                params_rev = ("v.a.xu@wustl.edu",password,f"ecNumber*{enz}", "ecNumber*", "substrates*", "commentarySubstrates*", "literatureSubstrates*", "organismSubstrates*", "products*", "commentaryProducts*", "literatureProducts*", "organismProducts*", "reversibility*")
                try:
                    res_rev = client.service.getSubstratesProducts(*params_rev)
                except Exception as e:
                    print('BRENDA request error:',e)
                    time.sleep(1/3)
                    res_rev = client.service.getSubstratesProducts(*params_rev)
                reversibility = [x['reversibility'] for x in res_nat_rev+res_rev]

                natRxnSubs = {}
                for idx,partner in enumerate(natRxnPartners):
                    if partner not in list(natRxnSubs.keys()):
                        natRxnSubs[partner+f'_{idx}'] = []
                    natRxnSubs[partner+f'_{idx}'].append(naturalSubstrate[idx])
                # pare down to unique values in values
                for key in list(natRxnSubs.keys()):
                    natRxnSubs[key] = list(set(natRxnSubs[key]))
                # pare down to unique values in keys
                natRxnSubs2 = {}
                for key in list(natRxnSubs.keys()):
                    key_2 = re.split('_',key)[0]
                    if key_2 not in list(natRxnSubs2.keys()):
                        natRxnSubs2[key_2] = []
                    natRxnSubs2[key_2] += [x.lower() for x in natRxnSubs[key] if ' H+' not in x and x != 'H+']
                # combining foward and reverse reaction partners
                natRxnSubs3 = {}
                for key_2 in list(natRxnSubs2.keys()):
                    split_list = [y.split(' + ') for y in [x.lower() for x in re.split(' = ',key_2)]]
                    replaced_list = [sorted([x.lower() for x in split_list[0] if ' H+' not in x and x != 'H+']),sorted([x.lower() for x in split_list[1] if ' H+' not in x and x != 'H+'])]
                    sorted_list = ['_'.join(x) for x in replaced_list]
                    key_3 = sorted_list[0]+'='+sorted_list[1]
                    if '?' not in key_3 and key_3 not in list(natRxnSubs3.keys()):
                        natRxnSubs3[key_3] = []
                    if '?' not in key_3:
                        natRxnSubs3[key_3] += natRxnSubs2[key_2]
                # pare down to unique values for each reaction partner
                for key_3 in list(natRxnSubs3.keys()):
                    natRxnSubs3[key_3] = list(set(natRxnSubs3[key_3]))
                # print('natRxnSubs3',natRxnSubs3)

                ### STARTHERE: need to account for specific reactions listed in BRENDA but not in KEGG
                ### need to figure out how to get the generic reactions to work
                natRxnSubs3_CIDs = {}
                for key_3 in list(natRxnSubs3.keys()):
                    split_list = [y.split('_') for y in [x for x in re.split('=',key_3)]]
                    sub_meta_names = [re.split(' ',x)[1] if re.split(' ',x)[0].isdigit() else x for x in [x.lower() for x in split_list[0]]]
                    prod_meta_names = [re.split(' ',x)[1] if re.split(' ',x)[0].isdigit() else x for x in [x.lower() for x in split_list[1]]]
                    for meta_name in sub_meta_names+prod_meta_names:
                        if not meta_name in list(self.__BRENDA_ligand_name_to_CID.keys()):
                            try:
                                print('meta_name',meta_name)
                                CID_req_raw = requests.get('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/'+meta_name+'/cids/JSON').text
                                parsed_CID_req = json.loads(CID_req_raw)
                                print('parsed_CID_req',parsed_CID_req)
                                CID_from_BRENDA = parsed_CID_req['IdentifierList']['CID']
                                self.__BRENDA_ligand_name_to_CID[meta_name] = CID_from_BRENDA
                            except Exception as e:
                                print('could not query correctly',e)

                    key_CIDs = ''
                    sub_key_CID_arr = []
                    prod_key_CID_arr = []
                    # print('meta_names',meta_names)
                    for idx,name in enumerate(sub_meta_names+prod_meta_names):
                        try:
                            if idx < len(sub_meta_names):
                                sub_key_CID_arr.append(self.__BRENDA_ligand_name_to_CID[name])
                            else:
                                prod_key_CID_arr.append(self.__BRENDA_ligand_name_to_CID[name])
                        except Exception as e:
                            print('meta not found in PubChem',e)

                    # remove H^+
                    if [1038] in sub_key_CID_arr:
                        sub_key_CID_arr.remove([1038])
                    if [1038] in prod_key_CID_arr:
                        prod_key_CID_arr.remove([1038])

                    CID_arr = sub_key_CID_arr+prod_key_CID_arr
                    CID_arr.sort()
                    for CID in CID_arr:
                        try:
                            if key_CIDs == '':
                                underscore = ''
                            else:
                                underscore = '_'
                            key_CIDs += underscore+str(CID)
                        except Exception as e:
                            print('meta not found in PubChem',e)

                    val_CIDs = []
                    val_meta_names = natRxnSubs3[key_3]
                    for val_meta_name in val_meta_names:
                        try:
                            val_CIDs += self.__BRENDA_ligand_name_to_CID[val_meta_name]
                        except Exception as e:
                            print('meta not found',e)
                    natRxnSubs3_CIDs[key_CIDs] = val_CIDs

                print('test6',natRxnSubs3)

                ### STARTHERE: figure out if all reactions listed in the two BRENDA queries (getNaturalSubstrates and getNaturalSubstratesProducts)
                # are being accounted for
                nat_rxn_partners_from_res_rev = {}
                for idx,res_item in enumerate(res_nat_rev+res_rev):
                    try:
                        res_item_py = zeep.helpers.serialize_object(res_item)
                        keys = list(res_item_py.keys())
                        if 'naturalSubstrates' in keys and 'naturalProducts' in keys:
                            res_item_2 = res_item_py['naturalSubstrates']
                            res_item_3 = res_item_py['naturalProducts']
                        elif 'substrates' in keys and 'products' in keys:
                            res_item_2 = res_item_py['substrates']
                            res_item_3 = res_item_py['products']
                        nat_subs = sorted([x.lower() for x in res_item_2.split(' + ') if ' H+' not in x and x != 'H+'])
                        nat_prods = sorted([x.lower() for x in res_item_3.split(' + ') if ' H+' not in x and x != 'H+'])
                        # print('nat_subs',nat_subs,nat_prods)
                        sub_meta_names = [re.split(' ',x)[1] if re.split(' ',x)[0].isdigit() else x for x in nat_subs]
                        prod_meta_names = [re.split(' ',x)[1] if re.split(' ',x)[0].isdigit() else x for x in nat_prods]

                        key = ''
                        sub_key_CID_arr = []
                        prod_key_CID_arr = []
                        for idx_2,meta_name in enumerate(sub_meta_names+prod_meta_names):
                            if idx_2 < len(sub_meta_names):
                                sub_key_CID_arr.append(self.__BRENDA_ligand_name_to_CID[meta_name])
                            else:
                                prod_key_CID_arr.append(self.__BRENDA_ligand_name_to_CID[meta_name])
                        CID_arr = sub_key_CID_arr+prod_key_CID_arr
                        CID_arr.sort()
                        for CID in CID_arr:
                            if key == '':
                                underscore = ''
                            else:
                                underscore = '_'
                            key += underscore+str(CID)
                        if '?' not in key:
                            # print('ttest5',key,idx,reversibility)
                            if key not in list(nat_rxn_partners_from_res_rev.keys()):
                                nat_rxn_partners_from_res_rev[key] = []
                            nat_rxn_partners_from_res_rev[key].append(reversibility[idx])
                            # print('nat_rxn_partners_5',nat_rxn_partners_from_res_rev)
                    except Exception as e:
                        exc_type, exc_obj, exc_tb = sys.exc_info()
                        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                        print(exc_type, fname, exc_tb.tb_lineno)
                        print('rxn_reversibility_err',e)

                for key in list(nat_rxn_partners_from_res_rev.keys()):
                    if 'r' in nat_rxn_partners_from_res_rev[key]:
                        nat_rxn_partners_from_res_rev[key] = 'r'
                    else:
                        nat_rxn_partners_from_res_rev[key] = 'ir'
                # print('test65',nat_rxn_partners_from_res_rev)

                # storing reversibility in each MPNG_Reaction
                keys = [x.replace(';','') for x in self.__enzymes[enz].reactions]
                rxn_reversible = {key:'both' for key in keys}
                # find MPNG_Reaction by checking if all metabolites are in the names of the metabolites in KEGG reaction
                for key_4 in self.__enzymes[enz].reactions:
                    key_4 = key_4.replace(';','')
                    print('key_4',key_4)
                    if 'G' not in key_4 and '(' not in key_4:
                        rxn_obj = self.__reactions[key_4]
                        sub_meta_entries = [x.id for x in list(rxn_obj.stoich.keys()) if rxn_obj.stoich[x] < 0 and x.id != 'C00080']
                        prod_meta_entries = [x.id for x in list(rxn_obj.stoich.keys()) if rxn_obj.stoich[x] > 0 and x.id != 'C00080']
                        print('prod_meta_entries',prod_meta_entries,sub_meta_entries)
                        try:
                            sub_meta_CIDs = [self.__MPNG_Metabolite_to_CID[x] for x in sub_meta_entries]
                            prod_meta_CIDs = [self.__MPNG_Metabolite_to_CID[x] for x in prod_meta_entries]
                            meta_CIDs = sub_meta_CIDs+prod_meta_CIDs
                            meta_CIDs.sort()
                            print('meta_CIDs',meta_CIDs)
                        except Exception as e:
                            print('meta name is generic',e)
                            continue
                        meta_CIDs_str = ''
                        for CID in sub_meta_CIDs+prod_meta_CIDs:
                            if meta_CIDs_str == '':
                                underscore = ''
                            else:
                                underscore = '_'
                            meta_CIDs_str += underscore+str(CID)

                        # print('test1241',natRxnSubs3_CIDs,meta_CIDs_str,meta_CIDs_str in list(natRxnSubs3_CIDs.keys()))
                        try:
                            print('nat',natRxnSubs3_CIDs)
                            natRxnSubs3_CIDs_arr = []
                            for CIDs_item in list(natRxnSubs3_CIDs.keys()):
                                CIDs_2_arr = []
                                for CID in CIDs_item.split('_'):
                                    CID_2 = CID.replace('[','').replace(']','')
                                    CIDs_3_arr = []
                                    for CID_3 in CID_2.split(', '):
                                        CIDs_3_arr.append(int(CID_3))
                                    CIDs_2_arr.append(CIDs_3_arr)
                                natRxnSubs3_CIDs_arr.append(CIDs_2_arr)
                            meta_CIDs_arr = []
                            for CID in meta_CIDs_str.split('_'):
                                CID_2 = CID.replace('[','').replace(']','')
                                CIDs_3_arr = []
                                for CID_3 in CID_2.split(', '):
                                    CIDs_3_arr.append(int(CID_3))
                                meta_CIDs_arr.append(CIDs_3_arr)

                            print('test2',natRxnSubs3_CIDs_arr,meta_CIDs_arr)
                            # print('test1',any([all([any([y in x for x in z]) for y in meta_CIDs_arr]) for z in natRxnSubs3_CIDs_arr]))
                            matching_rxn = False
                            rxn_item_idx = 0
                            for idx,CIDs_1 in enumerate(natRxnSubs3_CIDs_arr):
                                for CIDs_2 in meta_CIDs_arr:
                                    for x in CIDs_2:
                                        for y in CIDs_1:
                                            if x in y:
                                                CIDs_1.remove(y)
                                            if [] in natRxnSubs3_CIDs_arr:
                                                matching_rxn = True
                                                rxn_item_idx = idx
                                                # print('again',idx,CIDs_1,CIDs_2,x,y)
                                                break
                                            print('test2',natRxnSubs3_CIDs_arr,meta_CIDs_arr)
                                        else: continue
                                        break
                                    else: continue
                                    break
                                else: continue
                                break

                            if matching_rxn:
                            # if meta_CIDs_str in list(natRxnSubs3_CIDs.keys()):
                                print('meta_CIDs',meta_CIDs)
                                print('natRxnSubs3_CIDs',rxn_item_idx,natRxnSubs3_CIDs[list(natRxnSubs3_CIDs.keys())[rxn_item_idx]])
                                subs_rev = all([x in natRxnSubs3_CIDs[list(natRxnSubs3_CIDs.keys())[rxn_item_idx]] for x in meta_CIDs])
                                if list(natRxnSubs3_CIDs.keys())[rxn_item_idx] in list(nat_rxn_partners_from_res_rev.keys()):
                                    subs_prod_rev = nat_rxn_partners_from_res_rev[list(natRxnSubs3_CIDs.keys())[rxn_item_idx]] == 'r'
                                else:
                                    subs_prod_rev = False
                                tot_rev = subs_rev or subs_prod_rev
                            else:
                                tot_rev = True

                            if tot_rev:
                                rxn_reversible[key_4] = 'both'
                            else:
                                stoich = self.__reactions[key_4].stoich
                                subs_CIDs = []
                                prod_CIDs = []
                                for meta in list(stoich.keys()):
                                    try:
                                        # print('test6',self.__metabolites[meta.id].names)
                                        if stoich[meta] < 0:
                                            subs_CIDs += self.__MPNG_Metabolite_to_CID[meta.id]
                                        else:
                                            prod_CIDs += self.__MPNG_Metabolite_to_CID[meta.id]
                                    except Exception as e:
                                        print('metabolite not in system',e)

                                print('subs_CIDs',key_4,subs_CIDs,rxn_item_idx,list(natRxnSubs3_CIDs.keys())[rxn_item_idx],natRxnSubs3_CIDs[list(natRxnSubs3_CIDs.keys())[rxn_item_idx]],natRxnSubs3_CIDs)
                                if all([x in subs_CIDs for x in natRxnSubs3_CIDs[list(natRxnSubs3_CIDs.keys())[rxn_item_idx]]]):
                                    rxn_reversible[key_4] = 'forward'
                                else:
                                    rxn_reversible[key_4] = 'backward'
                        except Exception as e:
                            print('set rev err',e)

                print('rxn_reversible',rxn_reversible)

                for rxn_entry in self.__enzymes[enz].reactions:
                    rxn_entry = rxn_entry.replace(';','')
                    if 'G' not in key_4 and '(' not in key_4:
                        self.__reactions[rxn_entry].check_reversibility(enz,rxn_reversible)

            except Exception as e:
                print('enzyme '+enz+' did not work',e)

        print('updating reactions.json')
        with open('./src/stores/reactions.json', 'w', encoding='utf-8') as f:
            json.dump(list(map(lambda x: x.toJSON(),list(self.__reactions.values()))), f, ensure_ascii=False, indent=4)

    def minimize_side_reactions(self,graph:MetabolicPathwayNetworkGraph):

        return

    def calculate_kinetics(self,graph:MetabolicPathwayNetworkGraph):

        return

    def evaluate_reaction_kinetics(self,graph:MetabolicPathwayNetworkGraph):

        return

    def generate_whole_network(self,network_name:str,invalid_enz:list[str],invalid_rxn:list[str]) -> MetabolicPathwayNetworkGraph:
        # Task 1: construct metabolic network connections
        self.__graphs[network_name] = MetabolicPathwayNetworkGraph(network_name)

        # self.metabolites = self.__get_metabolites('all')

        for rxn in [x for x in self.__get_reactions('all') if (x.forward_valid or x.backward_valid) and x.entry not in invalid_rxn]:
            meta_ids = list(map(lambda y: y.id,list(rxn.stoich.keys())))
            enz_entries = list(rxn.enzyme_id.keys())
            matches = list(set(enz_entries) & set(invalid_enz))
            for match in matches:
                enz_entries.remove(match)
            if rxn.entry == 'R00753':
                print('worked',enz_entries)
            if all([z not in meta_ids for z in self.__excluded_metabolite_entries]) and len(enz_entries) > 0:
                try:
                    self.__graphs[network_name].add_reaction(rxn,[self.__get_metabolites(meta) for meta in list(map(lambda x: x.id,list(rxn.stoich.keys()))) if 'G' not in meta and '(' not in meta])
                except Exception as e:
                    print('meta not listed in KEGG',e)

        self.__graphs[network_name].generate_COBRA_model()

        return self.__graphs[network_name]

    def seek_optimal_network(self,
                             network_name:str,
                             objective_meta_entries:list[str],
                             substrate_metabolite_entries:list[str],
                             min_enzyme_ct:int,
                             max_enzyme_ct:int) -> None:
        # 1. Find optimal network via constrained mass balance
        self.__graphs[network_name].seek_optimal_network(objective_metabolite_entry=objective_meta_entries[0],
                                                     substrate_metabolite_entries=substrate_metabolite_entries,
                                                     min_enzyme_ct=min_enzyme_ct,
                                                     max_enzyme_ct=max_enzyme_ct)

        # Task 4: minimize side reactions
        self.minimize_side_reactions(self.__graphs[network_name])

        # Task 5: calculate kinetic parameters
        self.calculate_kinetics(self.__graphs[network_name])

        # Task 6: evaluate reaction kinetics
        self.evaluate_reaction_kinetics(self.__graphs[network_name])

        # Task 6: visualize results
        return self.__graphs[network_name]

    def visualize_graph(self,network_name:str) -> None:
        MPNG_net = self.__graphs[network_name]
        MPNG_net.vis_Network = Network()
        excluded_meta_entries = [x for x in self.__common_metabolite_entries+self.__small_gas_metas_entries if x != 'C00024']
        for idx,co_rxn in enumerate(MPNG_net.COBRA_model.reactions):
            flux_val = MPNG_net.mass_balance_sln.fluxes[co_rxn.id]
            # adding edges to slim graph
            if abs(round(flux_val)):
                if len([rxn for rxn in MPNG_net.get_reactions('all') if rxn.entry == co_rxn.id]) == 0:
                    continue
                rxn = [rxn for rxn in MPNG_net.get_reactions('all') if rxn.entry == co_rxn.id][0]
                if flux_val < 0:
                    for key in list(rxn.enzyme_id.keys()):
                        if rxn.enzyme_id[key] == 'backward' or rxn.enzyme_id[key] == 'both':
                            curr_enz_id = key
                            break
                else:
                    for key in list(rxn.enzyme_id.keys()):
                        if rxn.enzyme_id[key] == 'forward' or rxn.enzyme_id[key] == 'both':
                            curr_enz_id = key
                            break

                try:
                    curr_enz_name = curr_enz_id
                    # curr_enz_name = self.__get_enzymes(curr_enz_id).names[0]
                except Exception as e:
                    curr_enz_name = ''
                    print('no enzyme found for this reaction')

                if rxn.entry+'_'+curr_enz_name in MPNG_net.vis_Network.get_nodes():
                    MPNG_net.vis_Network.get_node(rxn.entry+'_'+curr_enz_name)['label'] = MPNG_net.vis_Network.get_node(rxn.entry+'_'+curr_enz_name)['label']+str('; enz_f: '+str(round(abs(flux_val))))
                else:
                    MPNG_net.vis_Network.add_node(rxn.entry+'_'+curr_enz_name,rxn.entry+'_'+curr_enz_name+'; enz_f: '+str(round(abs(flux_val))),shape='box')
                for meta in list(co_rxn.metabolites.keys()):
                    if meta.id in excluded_meta_entries:
                        MPNG_net.vis_Network.add_node(meta.id+'_'+str(idx),[metabolite.names[0] for metabolite in MPNG_net.get_metabolites('all') if metabolite.entry == meta.id][0],shape='image',image='https://rest.kegg.jp/get/'+meta.id+'/image')
                    else:
                        # print('graph_test',meta.id,[metabolite.names[0] for metabolite in MPNG_net.get_metabolites('all') if metabolite.entry == meta.id])
                        try:
                            MPNG_net.vis_Network.add_node(meta.id,[metabolite.names[0] for metabolite in MPNG_net.get_metabolites('all') if metabolite.entry == meta.id][0],shape='image',image='https://rest.kegg.jp/get/'+meta.id+'/image')
                        except Exception as e:
                            print('Glycan node in model, no name.',e)
                    try:
                        if co_rxn.metabolites[meta] < 0:
                            arrow_dir = 'to' if round(flux_val) > 0 else 'from'
                            if meta.id in excluded_meta_entries:
                                MPNG_net.vis_Network.add_edge(meta.id+'_'+str(idx),rxn.entry+'_'+curr_enz_name,label='rxn_f: '+str(int(abs(co_rxn.metabolites[meta]*flux_val))),arrows=arrow_dir,color='red')
                            else:
                                MPNG_net.vis_Network.add_edge(meta.id,rxn.entry+'_'+curr_enz_name,label='rxn_f: '+str(int(abs(co_rxn.metabolites[meta]*flux_val))),arrows=arrow_dir,width=3)
                        elif co_rxn.metabolites[meta] > 0:
                            arrow_dir = 'to' if round(flux_val) > 0 else 'from'
                            if meta.id in excluded_meta_entries:
                                MPNG_net.vis_Network.add_edge(rxn.entry+'_'+curr_enz_name,meta.id+'_'+str(idx),label='rxn_f: '+str(int(abs(co_rxn.metabolites[meta]*flux_val))),arrows=arrow_dir,color='red')
                            else:
                                MPNG_net.vis_Network.add_edge(rxn.entry+'_'+curr_enz_name,meta.id,label='rxn_f: '+str(int(abs(co_rxn.metabolites[meta]*flux_val))),arrows=arrow_dir,width=3)
                    except Exception as e:
                        print('Glycan edge in model, no name.',e)

        MPNG_net.vis_Network.layout = False
        MPNG_net.vis_Network.options.physics.enabled = True
        MPNG_net.vis_Network.show_buttons()
        MPNG_net.vis_Network.show(network_name+'_slim.html',local=True,notebook=False)
        return