import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from sympy import UnevaluatedExpr
import time
from cobra import Model, Reaction, Metabolite, Solution
from cobra.util.solver import add_absolute_expression, choose_solver, check_solver, get_solver_name
from optlang.symbolics import Zero
import networkx as nx
from pyvis.network import Network
import re
import equilibrator_api as eq

from src.classes.components.MPNG_Metabolite import MPNG_Metabolite
from src.classes.components.MPNG_Reaction import MPNG_Reaction
from src.classes.components.MPNG_Enzyme import MPNG_Enzyme

class MetabolicPathwayNetworkGraph:
    def __init__(self,name:str) -> None:
        self.__name = name
        self.__metabolites: dict[str,MPNG_Metabolite] = {}
        self.__reactions: dict[str,MPNG_Reaction] = {}
        self.__enzymes: dict[str,MPNG_Enzyme] = {}

        self.__COBRA_model: Model = Model(name+'_cobra')
        self.__mass_balance_sln: Solution = None
        self.__rxn_dGr: dict = {}
        self.__all_possible_reactions: list[str] = []
        self.__vis_Graph: nx.Graph = nx.DiGraph()
        self.__path_Graph: nx.DiGraph = nx.DiGraph()

        self.__temperature = 303.15 # K
        self.__pH = 7
        self.E_carriers = {
            'ATP': 0, 'ADP': 0,
            'NADH': 0, 'NAD': 0,
            'NADPH': 0, 'NADP': 0,
            'FADH2': 0, 'FAD': 0
        }
        self.__cc = eq.ComponentContribution()

        self.__small_liquid_meta_entries = ['C00001']
        self.__small_gas_meta_entries = ['C00007','C00011']
        self.__common_metabolite_entries = ['C00138','C00139','C00080','C00024','C00125','C00126']
        for x in range(14):
            zeros = '0'*(5-len(str(x+1)))
            self.__common_metabolite_entries.append('C'+zeros+str(x+1))
        # self.__common_metabolites = self.get_metabolites(self.__common_metabolite_entries)

    @property
    def name(self) -> str:
        return self.__name

    @property
    def COBRA_model(self) -> Model:
        return self.__COBRA_model

    @COBRA_model.setter
    def COBRA_model(self,new_model:Model) -> None:
        self.__COBRA_model = new_model

    @property
    def all_possible_reactions(self) -> Model:
        return self.__all_possible_reactions

    @property
    def vis_Network(self) -> nx.DiGraph:
        return self.__vis_Graph

    @vis_Network.setter
    def vis_Network(self,new_graph:nx.DiGraph) -> None:
        self.__vis_Graph = new_graph

    @property
    def path_Graph(self) -> nx.DiGraph:
        return self.__path_Graph

    @path_Graph.setter
    def path_Graph(self,new_graph:nx.DiGraph) -> None:
        self.__path_Graph = new_graph

    @property
    def root_metabolite(self) -> MPNG_Metabolite:
        return self.__root_metabolite

    @property
    def leaf_metabolites(self) -> list[MPNG_Metabolite]:
        return self.__leaf_metabolites
    
    @leaf_metabolites.setter
    def leaf_metabolites(self,leaves:list[MPNG_Metabolite]) -> None:
        self.__leaf_metabolites = leaves

    @property
    def mass_balance_sln(self) -> Solution:
        return self.__mass_balance_sln

    @mass_balance_sln.setter
    def mass_balance_sln(self,new_sln:Solution) -> None:
        self.__mass_balance_sln = new_sln

    @property
    def metabolites(self) -> dict[str,MPNG_Metabolite]:
        return self.__metabolites

    @metabolites.setter
    def metabolites(self,new_metabolites:MPNG_Metabolite|list[MPNG_Metabolite]|list) -> None:
        if type(new_metabolites) is MPNG_Metabolite:
            self.__metabolites[str(new_metabolites.entry)] = new_metabolites
        elif type(new_metabolites) is list:
            for m in new_metabolites:
                self.__metabolites[str(m.entry)] = m
        elif type(new_metabolites) is list and new_metabolites == []:
            self.__metabolites = {}

    def remove_meta_by_entry(self,entries:list[str]) -> None:
        for x in entries:
            del self.__metabolites[x]

    def get_metabolites(self,entries:str|list) -> MPNG_Metabolite | list[MPNG_Metabolite]:
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

    @property
    def reactions(self) -> dict[str,MPNG_Reaction]:
        return self.__reactions

    @reactions.setter
    def reactions(self,new_reactions:MPNG_Reaction|list) -> None:
        if type(new_reactions) is MPNG_Reaction:
            self.__reactions[str(new_reactions.entry)] = new_reactions
        elif type(new_reactions) is list:
            # self.__reactions = {}
            for r in new_reactions:
                self.__reactions[str(r.entry)] = r
        elif type(new_reactions) is list and new_reactions == []:
            self.__reactions = {}

    def remove_rxns_by_entry(self,entries:list[str]) -> None:
        for x in entries:
            del self.__reactions[x]

    def get_reactions(self,entries:str|list) -> MPNG_Reaction | list[MPNG_Reaction]:
        if type(entries) == str:
            if entries == 'all':
                return list(self.__reactions.values())
            else:
                return self.__reactions[entries]
        elif type(entries) == list:
            rxns = []
            # print('testest',entries)
            for x in entries:
                # print('testest2',self.__reactions[x])
                rxns.append(self.__reactions[x])
            return rxns
        else:
            return []

    @property
    def enzymes(self) -> dict[str,MPNG_Enzyme]:
        return self.__enzymes

    @enzymes.setter
    def enzymes(self,new_enzymes:MPNG_Enzyme|list[MPNG_Enzyme]|list) -> None:
        if type(new_enzymes) is MPNG_Enzyme:
            self.__enzymes[str(new_enzymes.entry)] = new_enzymes
        elif type(new_enzymes) is list:
            self.__enzymes = {}
            for r in new_enzymes:
                self.__enzymes[str(r.entry)] = r
        elif type(new_enzymes) is list and new_enzymes == []:
            self.__enzymes = {}

    def get_enzymes(self,entries:str|list) -> MPNG_Enzyme | list[MPNG_Enzyme]:
        if type(entries) == str:
            if entries == 'all':
                return list(self.__enzymes.values())
            else:
                return self.__enzymes[entries]
        elif type(entries) == list:
            enz = []
            for x in entries:
                enz.append(self.__enzymes[x])
        else:
            return []

    @property
    def temperature(self) -> float:
        return self.__temperature

    @temperature.setter
    def temperature(self,new_temp:float) -> None:
        self.__temperature = new_temp

    @property
    def pH(self) -> float:
        return self.__pH

    @pH.setter
    def pH(self,new_pH:float) -> None:
        self.__pH = new_pH

    def __update_COBRA_model(self,new_reactions:list[MPNG_Reaction]) -> None:
        # COBRA automatically checks if reaction already exists (ignored if it does)
        new_rxns = []
        for rxn in new_reactions:
            ub = 1E6
            lb = -1E6
            if not rxn.forward_valid:
                ub = 0
            if not rxn.backward_valid:
                lb = 0
            reaction: Reaction = Reaction(rxn.entry,lower_bound=lb,upper_bound=ub)
            reaction.add_metabolites(rxn.stoich)
            new_rxns.append(reaction)
        print('adding COBRA reactions...')
        self.__COBRA_model.add_reactions(new_rxns)
        print('all COBRA reactions added')

    def add_reaction(self,new_reaction:MPNG_Reaction,new_metabolites:list[MPNG_Metabolite]) -> None:
        # add MPNG_Reaction to MPNG
        self.reactions = new_reaction
        # print('add_reaction',new_reaction.enzyme_id)
        # rxn_number = new_reaction.enzyme_id[0]+':'+new_reaction.entry+'_'+str(len(self.__reactions.keys()))
        # add to NX Graph
        self.__vis_Graph.add_node(node_for_adding=new_reaction.entry,id=len(self.__vis_Graph.nodes))
        self.__path_Graph.add_node(node_for_adding=new_reaction.entry+'_f',id=len(self.__path_Graph.nodes))
        self.__path_Graph.add_node(node_for_adding=new_reaction.entry+'_r',id=len(self.__path_Graph.nodes))

        stoich: dict[Metabolite,int] = new_reaction.stoich
        key_entries = list(map(lambda x: x.id,list(stoich.keys())))
        new_metabolite_entries = list(map(lambda x: x.entry,new_metabolites))
        for idx,m in enumerate(new_metabolite_entries):
            if m not in self.__vis_Graph.nodes and 'G' not in m and '(' not in m:
                self.__vis_Graph.add_node(node_for_adding=m)
                self.metabolites = new_metabolites[idx]
                if stoich[list(stoich.keys())[key_entries.index(m)]] < 0:
                    self.__vis_Graph.add_edge(m,new_reaction.entry,
                        stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    )
                elif stoich[list(stoich.keys())[key_entries.index(m)]] > 0:
                    self.__vis_Graph.add_edge(new_reaction.entry,m,
                        stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    )
                    # self.__path_Graph.add_edge(new_reaction.entry+'_f',m,
                    #     stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    # )
                    # self.__path_Graph.add_edge(m,new_reaction.entry+'_r',
                    #     stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    # )
            if 'G' not in m and '(' not in m and m not in [x for x in self.__common_metabolite_entries if x != 'C00011']:
                self.__path_Graph.add_node(node_for_adding=m)
                if stoich[list(stoich.keys())[key_entries.index(m)]] < 0:
                    self.__path_Graph.add_edge(m,new_reaction.entry+'_f',
                        stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    )
                    self.__path_Graph.add_edge(new_reaction.entry+'_r',m,
                        stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    )
                elif stoich[list(stoich.keys())[key_entries.index(m)]] > 0:
                    self.__path_Graph.add_edge(new_reaction.entry+'_f',m,
                        stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    )
                    self.__path_Graph.add_edge(m,new_reaction.entry+'_r',
                        stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                    )

    def add_reaction_slim(self,new_reaction:MPNG_Reaction,all_metas:dict[str,MPNG_Metabolite]) -> None:
        # add to NX Graph
        self.__vis_Graph.add_node(node_for_adding=new_reaction.entry)
        self.__vis_Graph.add_node(node_for_adding=new_reaction.entry)

        # add to MPNG
        self.reactions = new_reaction

        stoich: dict[Metabolite,int] = new_reaction.stoich
        key_entries = list(map(lambda x: x.id,list(stoich.keys())))
        for m in key_entries:
            if m not in self.__vis_Graph.nodes:
                self.__vis_Graph.add_node(node_for_adding=m)
            if 'G' not in m and '(' not in m:
                self.metabolites = all_metas[m]
            if stoich[list(stoich.keys())[key_entries.index(m)]] < 0:
                self.__vis_Graph.add_edge(m,new_reaction.entry,
                    stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                )
            elif stoich[list(stoich.keys())[key_entries.index(m)]] > 0:
                self.__vis_Graph.add_edge(new_reaction.entry,m,
                    stoich=stoich[list(stoich.keys())[key_entries.index(m)]]
                )

    def generate_COBRA_model(self) -> None:
        print('generating COBRA model...')
        self.__update_COBRA_model(self.get_reactions('all'))
        print('generating COBRA boundaries...')
        boundary_rxns = [rxn for rxn in self.__COBRA_model.reactions if 'SK_' in rxn.id]
        self.__COBRA_model.remove_reactions(boundary_rxns)
        for meta in self.__COBRA_model.metabolites:
            try:
                self.__COBRA_model.add_boundary(meta,type='sink',lb=0,ub=0)
            except Exception as e:
                pass
        print('COBRA model ready')

    # what are the balancing metabolites that we can introduce?
    # N: helps produce C+
    #   N2, NH4
    # C: helps produce C-
    #   CO2
    #   Fatty acids?
    # O: helps produce C+
    #   H2O

    # is there a way to exclude H2O from the mass balance?
    def __calculate_redox_balance(self,objective_meta_entry:str,substrate_meta_entries:list[str]):
        balancing_meta_entries = ['C00001','C00014','C00011','C00697']
        redox_vals = {}
        element_cts = []
        for entry in [objective_meta_entry]+substrate_meta_entries+balancing_meta_entries:
            meta = self.get_metabolites(entry)
            formula_str = meta.formula

            element_ct = {'C_ct':0,'H_ct':0,'O_ct':0,'N_ct':0}
            letters = re.findall(r'[A-Z]',formula_str)
            nums = re.findall(r'\d+', formula_str)
            if len(letters) > len(nums):
                for idx,x in enumerate(letters):
                    letter_idx = formula_str.find(x)
                    if letter_idx+1 < len(formula_str) and not formula_str[letter_idx+1].isnumeric():
                        nums.insert(idx,'1')
                if not formula_str[-1].isnumeric():
                    nums.append('1')
            for idx,x in enumerate(letters):
                element_ct[x+'_ct'] = int(nums[idx])
            redox_vals[entry] = 2*element_ct['O_ct'] - element_ct['H_ct'] + 3*element_ct['N_ct']
            element_cts.append(element_ct)

        if sum([x['N_ct'] for x in element_cts[0:len([objective_meta_entry]+substrate_meta_entries)]]) == 0:
            balancing_meta_entries = ['C00001','C00011']

            def func(x):
                res = [0]*6
                # mass and redox balance the reaction
                for idx,entry in enumerate([objective_meta_entry]+substrate_meta_entries+balancing_meta_entries):
                    el_ct = element_cts[idx]

                    # C mass balance
                    res[0] += el_ct['C_ct']*x[idx]
                    # H mass balance
                    res[1] += el_ct['H_ct']*x[idx]
                    # O mass balance
                    res[2] += el_ct['O_ct']*x[idx]
                    # redox balance
                    res[3] += redox_vals[entry]*x[idx]

                res[4] = x[0] - 1000
                res[5] = x[-2] - 500

                return np.array(res)

            res = least_squares(func,[1000]+[-1000]*len(substrate_meta_entries)+[0]*len(balancing_meta_entries),
                                bounds=([0]+[-np.inf]*len(substrate_meta_entries)+[-np.inf]+[0],[np.inf]+[0]*len(substrate_meta_entries)+[np.inf]*len(balancing_meta_entries)))
        else:
            # balancing_meta_entries = ['C00014','C00011','C00697']

            def func(x):
                res = [0]*6
                # mass and redox balance the reaction
                for idx,entry in enumerate([objective_meta_entry]+substrate_meta_entries+balancing_meta_entries):
                    el_ct = element_cts[idx]

                    # C mass balance
                    res[0] += el_ct['C_ct']*x[idx]
                    # H mass balance
                    res[1] += el_ct['H_ct']*x[idx]
                    # O mass balance
                    res[2] += el_ct['O_ct']*x[idx]
                    # N mass balance
                    res[3] += el_ct['N_ct']*x[idx]
                    # redox balance
                    res[4] += redox_vals[entry]*x[idx]

                res[5] = x[0] - 1000

                return np.array(res)

            res = least_squares(func,[1000]+[-1000]*len(substrate_meta_entries)+[0]*len(balancing_meta_entries),
                                bounds=([0]+[-np.inf]*len(substrate_meta_entries)+[-np.inf]*2+[0]*2,[np.inf]+[0]*len(substrate_meta_entries)+[np.inf]*len(balancing_meta_entries)))

        print(res)
        overall_stoich = {key: res['x'][idx] for idx,key in enumerate([objective_meta_entry]+substrate_meta_entries+balancing_meta_entries)}

        return overall_stoich

    def __find_shortest_simple_paths(self,
                                     objective_meta_entry:str,
                                     substrate_meta_entries:list[str],
                                     min_enzyme_ct:int,
                                     max_enzyme_ct:int):
        shortest_paths = [x for x in nx.all_simple_paths(self.__path_Graph,source=objective_meta_entry,target=substrate_meta_entries,cutoff=max_enzyme_ct*2+1)]
        shortest_paths = [x for x in shortest_paths if len(x) >= min_enzyme_ct*2+1]
        return shortest_paths

    def __balance_candidate_networks(self,
                                     candidate_paths:list[list[str]],
                                     objective_meta_entry:str,
                                     substrate_meta_entries:list[str],
                                     overall_stoich:list[str]):
        # 1. Try to generate any solution from COBRA model while minimizing total flux - can modify this to use different constraint to obtain optimal solution
        print('number of candidate paths: ',len(candidate_paths),candidate_paths)
        start_1 = time.time()
        with self.__COBRA_model as mdl:
            for key in overall_stoich.keys():
                overall_stoich[key] = 10*(round(overall_stoich[key]/10))

            print(pd.DataFrame(list(overall_stoich.items()), columns=['meta entry','flux']))

            for meta_entry in list(overall_stoich.keys()):
                if overall_stoich[meta_entry] > 0:
                    mdl.reactions.get_by_id('SK_'+meta_entry).upper_bound = overall_stoich[meta_entry]
                    mdl.reactions.get_by_id('SK_'+meta_entry).lower_bound = overall_stoich[meta_entry]
                else:
                    mdl.reactions.get_by_id('SK_'+meta_entry).lower_bound = overall_stoich[meta_entry]
                    mdl.reactions.get_by_id('SK_'+meta_entry).upper_bound = overall_stoich[meta_entry]

            # shutting down reactions that utilize small_gas_metas as substrates - move this to initialization too?
            for rxn in [x for x in mdl.reactions if x not in mdl.boundary]:
                stoich = rxn.metabolites
                subs = []
                prods = []
                for key in list(stoich.keys()):
                    if stoich[key] > 0:
                        prods.append(key.id)
                    elif stoich[key] < 0:
                        subs.append(key.id)
                if any([x in prods for x in self.__small_gas_meta_entries]):
                    mdl.reactions.get_by_id(rxn.id).lower_bound = 0
                elif any([x in subs for x in self.__small_gas_meta_entries]):
                    mdl.reactions.get_by_id(rxn.id).upper_bound = 0

            # limiting exchanges to only those for objective and substrate metabolites and small gas metas
            for rxn in [x for x in mdl.boundary if re.split('_',x.id)[1] in self.__small_gas_meta_entries+self.__small_liquid_meta_entries]:
                if re.split('_',rxn.id)[1] in self.__small_gas_meta_entries:
                    mdl.reactions.get_by_id(rxn.id).lower_bound = 0
                else:
                    mdl.reactions.get_by_id(rxn.id).lower_bound = -1E6
                # only prevents intake of non-substrate metabolites, doesn't prevent production of other metabolites
                mdl.reactions.get_by_id(rxn.id).upper_bound = 1E6

            obj_dict = {}
            # 2.1 Create optimization objective for reactions on path
            path_opt_vars = []
            rxns_lvl_2 = []
            for idx,path in enumerate(candidate_paths):
                print('Creating optimization objective for path #',idx+1)
                for node in path:
                    if 'R' in node:
                        [rxn_entry,path_direction] = re.split('_',node)
                        rxn: Reaction = mdl.reactions.get_by_id(rxn_entry)
                        # here, the optimization stoichiometry is flipped because the pathway seeking starts from target  
                        # and seeks towards substrates
                        ### STARTHERE: try 0
                        if path_direction == 'r':
                            obj_dict[rxn.forward_variable] = 1
                            obj_dict[rxn.reverse_variable] = -1
                        else:
                            obj_dict[rxn.forward_variable] = -1
                            obj_dict[rxn.reverse_variable] = 1
                        path_opt_vars.append(rxn.forward_variable)
                        path_opt_vars.append(rxn.reverse_variable)
                        rxns_lvl_2.append(rxn)

            # 2.4 Create optimization objective for number of total sum of flux absolute values
            tot_flux_vars = []
            tot_flux_consts = []
            for rxn in [x for x in mdl.reactions if x.id not in list(map(lambda y: 'SK_'+y,self.__small_gas_meta_entries))]:
                new_var = mdl.problem.Variable('rxn_var_'+rxn.id)
                new_constraint = mdl.problem.Constraint(rxn.forward_variable + rxn.reverse_variable - new_var,
                                        name='rxn_constraint_'+rxn.id,
                                        ub=0,
                                        lb=0)
                tot_flux_vars.append(new_var)
                tot_flux_consts.append(new_constraint)
                mdl.add_cons_vars([new_var,new_constraint])
                obj_dict[new_var] = -1

            # 2.5 Preventing reactions that utilize target metabolite
            target_meta = mdl.metabolites.get_by_id(objective_meta_entry)
            for rxn in [x for x in mdl.reactions if x.id != 'SK_'+objective_meta_entry and x not in mdl.boundary]:
                stoich = rxn.metabolites
                subs = []
                prods = []
                for key in list(stoich.keys()):
                    if stoich[key] > 0:
                        prods.append(key)
                    elif stoich[key] < 0:
                        subs.append(key)
                if target_meta in prods:
                    mdl.reactions.get_by_id(rxn.id).lower_bound = 0
                elif target_meta in subs:
                    mdl.reactions.get_by_id(rxn.id).upper_bound = 0

            # 3.2 Generating optimization target based on reaction dG
            # and removing minimization of total flux
            dGr_dict = {}
            for rxn in rxns_lvl_2:
                dGr_dict[rxn] = -1

            for rxn in rxns_lvl_2:
                # method 1: minimize sum of abs dGr in network
                # method 2: if target metabolite is lower dGf than substrates, then try to minimize differences of dGf of intermediates
                # relative to substrates - do this by minimizing the value of dGr immediate to substrates and the sum of abs of the dGr of the whole network
                rxn_meta_keys = list(rxn.metabolites.keys())
                compound_ids = [meta.id for meta in rxn_meta_keys]
                rxn_dict = {}
                print('rxn_entry',rxn.id)
                for meta in rxn_meta_keys:
                    rxn_dict[self.__cc.get_compound(f"kegg:{meta.id}")] = rxn.metabolites[rxn_meta_keys[compound_ids.index(meta.id)]]
                    # print('rxn_2',rxn_meta_keys,compound_ids.index(meta.id),self.__cc.get_compound(f"kegg:{meta.id}"),rxn.metabolites[rxn_meta_keys[compound_ids.index(meta.id)]])
                try:
                    eqi_rxn = eq.Reaction(rxn_dict)
                    dGr = abs(self.__cc.dg_prime(eqi_rxn).value.m_as("kJ/mol"))
                    dGr_dict[rxn] = dGr
                    print('dGr',dGr)
                except Exception as e:
                    print('equilibrator warning: ',e)
                    print('abs dGr is assumed negative here')

            self.__rxn_dGr = dGr_dict

            if len(rxns_lvl_2) > 0:
                dGr_max = max(list(dGr_dict.values()))
                for rxn in list(dGr_dict.keys()):
                    if dGr_dict[rxn] >= 0:
                        dGr_val = dGr_dict[rxn]/dGr_max
                    else:
                        dGr_val = 0
                    print('test',dGr_val)
                    new_var = mdl.problem.Variable('abs_dGr_'+rxn.id)
                    new_constraint = mdl.problem.Constraint(rxn.forward_variable*dGr_val + rxn.reverse_variable*dGr_val - new_var,
                                            name='abs_dGr_constraint_'+rxn.id,
                                            ub=0,
                                            lb=0)
                    mdl.add_cons_vars([new_var,new_constraint])
                    obj_dict[new_var] = -1

            # 2.7 Print COBRA model metrics and solve level 1
            start_2 = time.time()
            print('obj_dict lvl_1',obj_dict)
            mdl.objective = mdl.problem.Objective(Zero, sloppy=True, direction="max")
            mdl.solver.objective.set_linear_coefficients(obj_dict)
            lvl_1_res = self.mass_balance_sln = mdl.optimize()
            print('objective_value lvl_1',lvl_1_res.objective_value)
            start_3 = time.time()

            print('start_1',start_2-start_1)
            print('start_2',start_3-start_2)

            fluxes_lvl_1 = self.mass_balance_sln.fluxes
            int_fluxes_lvl_1 = fluxes_lvl_1.loc[[x for x in fluxes_lvl_1.index if 'SK_' not in x]]
            bnd_fluxes_lvl_1 = fluxes_lvl_1.loc[[x for x in fluxes_lvl_1.index if 'SK_' in x]]
            rel_int_fluxes_lvl_1 = int_fluxes_lvl_1.loc[[x for x in int_fluxes_lvl_1.index if round(int_fluxes_lvl_1[x]) != 0]]
            rel_bnd_fluxes_lvl_1 = bnd_fluxes_lvl_1.loc[[x for x in bnd_fluxes_lvl_1.index if round(bnd_fluxes_lvl_1[x]) != 0]]
            print('number of internal rxns used (lvl_1): ',(round(int_fluxes_lvl_1)!=0).sum())
            print('number of boundary rxns used (lvl_1): ',(round(bnd_fluxes_lvl_1)!=0).sum())
            print('int_fluxes (lvl_1): ',rel_int_fluxes_lvl_1.to_markdown())
            print('bnd_fluxes (lvl_1): ',rel_bnd_fluxes_lvl_1.to_markdown())

            print('starting optimization lvl_3...')
            rxns_lvl_3 = [x for x in mdl.reactions if x.id in list(rel_int_fluxes_lvl_1.index)]
            # for rxn in [x for x in mdl.reactions if x not in rxns_lvl_3 and x not in mdl.boundary]:
            #     mdl.reactions.get_by_id(rxn.id).lower_bound = 0
            #     mdl.reactions.get_by_id(rxn.id).upper_bound = 0
            print('finished restricting to relevant reactions')

            # 4.1 choose the production pathway of target meta that generates the highest yield for each metabolite
            max_target_meta_prod_rxn = {}
            for rxn in [x for x in rxns_lvl_3 if 'SK_' not in x.id]:
                if objective_meta_entry in list(map(lambda y: y.id,list(rxn.metabolites.keys()))):
                    max_target_meta_prod_rxn[rxn.id] = abs(rel_int_fluxes_lvl_1.loc[rxn.id])
                    print('max',max_target_meta_prod_rxn)

            max_prod_rxn_ids = list(max_target_meta_prod_rxn.keys())
            print('max_prod_rxn_ids',max_prod_rxn_ids)
            if len(max_prod_rxn_ids) > 1:
                max_flux_rxn_id = max(max_target_meta_prod_rxn,key=max_target_meta_prod_rxn.get)
                print('test3',max_flux_rxn_id)
                non_maxed_target_fluxes = [x for x in max_prod_rxn_ids if x != max_flux_rxn_id]
                print('max_prod_rxn_2',non_maxed_target_fluxes)

                for rxn_id in non_maxed_target_fluxes:
                    mdl.reactions.get_by_id(rxn_id).upper_bound = 0
                    mdl.reactions.get_by_id(rxn_id).lower_bound = 0

                # 4.3 Print COBRA model metrics and solve level 3
                print('obj_dict lvl_3',obj_dict)
                mdl.objective = mdl.problem.Objective(Zero, sloppy=True, direction="max")
                mdl.solver.objective.set_linear_coefficients(obj_dict)
                lvl_3_res = self.mass_balance_sln = mdl.optimize()
                print('objective_value lvl_3',lvl_3_res.objective_value)

                fluxes_lvl_3 = self.mass_balance_sln.fluxes
                int_fluxes_lvl_3 = fluxes_lvl_3.loc[[x for x in fluxes_lvl_3.index if 'SK_' not in x]]
                bnd_fluxes_lvl_3 = fluxes_lvl_3.loc[[x for x in fluxes_lvl_3.index if 'SK_' in x]]
                rel_int_fluxes_lvl_3 = int_fluxes_lvl_3.loc[[x for x in int_fluxes_lvl_3.index if round(int_fluxes_lvl_3[x]) != 0]]
                rel_bnd_fluxes_lvl_3 = bnd_fluxes_lvl_3.loc[[x for x in bnd_fluxes_lvl_3.index if round(bnd_fluxes_lvl_3[x]) != 0]]
                print('number of internal rxns used (lvl_3): ',(round(int_fluxes_lvl_3)!=0).sum())
                print('number of boundary rxns used (lvl_3): ',(round(bnd_fluxes_lvl_3)!=0).sum())
                print('int_fluxes (lvl_3): ',rel_int_fluxes_lvl_3.to_markdown())
                print('bnd_fluxes (lvl_3): ',rel_bnd_fluxes_lvl_3.to_markdown())
        return

    def __calculate_equilibrium(self,objective_meta_entry:str,
                                          substrate_meta_entries:list[str],
                                          overall_stoich:pd.DataFrame):
        
        return

    def seek_optimal_network(self,
                             objective_metabolite_entry:str,
                             substrate_metabolite_entries:list[str],
                             min_enzyme_ct:int,
                             max_enzyme_ct:int) -> Solution:
        has_paths = []
        for meta in [x for x in substrate_metabolite_entries if x not in self.__common_metabolite_entries]:
            has_paths.append(nx.has_path(self.__path_Graph,meta,objective_metabolite_entry))
        if not any(has_paths):
            print('No simple path found between substrates and objective metabolite, exiting function.')
            return

        overall_stoich = self.__calculate_redox_balance(objective_meta_entry=objective_metabolite_entry,
                                                        substrate_meta_entries=substrate_metabolite_entries)

        # find potential paths between substrate(s) and objective metabolite
        candidate_paths = self.__find_shortest_simple_paths(objective_meta_entry=objective_metabolite_entry,
                                                            substrate_meta_entries=substrate_metabolite_entries,
                                                            min_enzyme_ct=min_enzyme_ct,
                                                            max_enzyme_ct=max_enzyme_ct)

        # balance all candidate paths
        self.__balance_candidate_networks(candidate_paths=[],
                                          objective_meta_entry=objective_metabolite_entry,
                                          substrate_meta_entries=substrate_metabolite_entries,
                                          overall_stoich=overall_stoich)

        # calculate equilibrium of all metabolites in system
        self.__calculate_equilibrium(objective_meta_entry=objective_metabolite_entry,
                                          substrate_meta_entries=substrate_metabolite_entries,
                                          overall_stoich=overall_stoich)

        return self.__COBRA_model