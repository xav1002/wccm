import numpy as np
import pandas as pd
from scipy import optimize as opt
from sympy import UnevaluatedExpr
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
        self.__temp_leaves = []
        self.__metabolites: dict[str,MPNG_Metabolite] = {}
        self.__reactions: dict[str,MPNG_Reaction] = {}
        self.__enzymes: dict[str,MPNG_Enzyme] = {}

        self.__COBRA_model: Model = Model(name+'_cobra')
        self.__mass_balance_sln: Solution = None
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
            # self.__metabolites = {}
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
            reaction: Reaction = Reaction(rxn.entry,lower_bound=None,upper_bound=None)
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

    def calc_dGr_explore_weight(self):
        return

    def calc_dHr_explore_weight(self):
        return

    def calc_rxn_redox_explore_weight(self):
        return

    def find_equilibrium_concentrations(self) -> None:
        return

    def generate_COBRA_model(self) -> None:
        print('generating COBRA model...')
        self.__update_COBRA_model(self.get_reactions('all'))
        print('generating COBRA boundaries...')
        boundary_rxns = [rxn for rxn in self.__COBRA_model.reactions if 'SK_' in rxn.id]
        self.__COBRA_model.remove_reactions(boundary_rxns)
        for meta in self.__COBRA_model.metabolites:
            try:
                self.__COBRA_model.add_boundary(meta,type='sink',lb=-1000,ub=1000)
            except Exception as e:
                pass
        print('COBRA model ready')

    def __find_shortest_simple_paths(self,
                                     objective_meta_entry:str,
                                     substrate_meta_entries:list[str],
                                     min_enzyme_ct:int,
                                     max_enzyme_ct:int):
        shortest_paths = [x for x in nx.all_simple_paths(self.__path_Graph,source=objective_meta_entry,target=substrate_meta_entries,cutoff=max_enzyme_ct*2+1)]
        shortest_paths = [x for x in shortest_paths if len(x) >= min_enzyme_ct*2+1]
        return shortest_paths

    def __balance_candidate_networks(self,candidate_paths:list[list[str]],objective_meta_entry:str,substrate_meta_entries:list[str]):
        # 1. Try to generate any solution from COBRA model while minimizing total flux - can modify this to use different constraint to obtain optimal solution
        bnd_meta_entries = substrate_meta_entries+[objective_meta_entry]
        print('number of candidate paths: ',len(candidate_paths),candidate_paths)
        with self.__COBRA_model as mdl:
            for rxn in [x for x in mdl.boundary if re.split('_',x.id)[1] in self.__small_liquid_meta_entries]:
                mdl.reactions.get_by_id(rxn.id).lower_bound = -10000
                mdl.reactions.get_by_id(rxn.id).upper_bound = 10000

            # Shutting down reverse reactions that are infeasible according to BRENDA
            for rxn in [x for x in mdl.reactions if x not in mdl.boundary]:
                # print('test',rxn.id,self.__reactions[rxn.id].forward_valid,self.__reactions[rxn.id].backward_valid)
                if self.__reactions[rxn.id].forward_valid:
                    mdl.reactions.get_by_id(rxn.id).upper_bound = 1000
                else: mdl.reactions.get_by_id(rxn.id).upper_bound = 0
                if self.__reactions[rxn.id].backward_valid:
                    mdl.reactions.get_by_id(rxn.id).lower_bound = -1000
                else: mdl.reactions.get_by_id(rxn.id).lower_bound = 0

            # shutting down reactions that utilize small_gas_metas as substrates
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
            for rxn in [x for x in mdl.boundary if re.split('_',x.id)[1] not in bnd_meta_entries and re.split('_',x.id)[1] not in self.__small_gas_meta_entries+self.__small_liquid_meta_entries]:
                mdl.reactions.get_by_id(rxn.id).lower_bound = 0
                # only prevents intake of non-substrate metabolites, doesn't prevent production of other metabolites
                mdl.reactions.get_by_id(rxn.id).upper_bound = 0

            for rxn in [x for x in mdl.boundary if re.split('_',x.id)[1] in self.__small_gas_meta_entries]:
                mdl.reactions.get_by_id(rxn.id).lower_bound = 0
                mdl.reactions.get_by_id(rxn.id).upper_bound = 10000

            obj_dict = {}
            # 2.1 Create optimization objective for reactions on path
            path_opt_vars = []
            for idx,path in enumerate(candidate_paths):
                print('Creating optimization objective for path #',idx+1)
                for node in path:
                    if 'R' in node:
                        [rxn_entry,path_direction] = re.split('_',node)
                        rxn: Reaction = mdl.reactions.get_by_id(rxn_entry)
                        # here, the optimization stoichiometry is flipped because the pathway seeking starts from target  
                        # and seeks towards substrates
                        if path_direction == 'r':
                            obj_dict[rxn.forward_variable] = 2
                            obj_dict[rxn.reverse_variable] = -2
                        else:
                            obj_dict[rxn.forward_variable] = -2
                            obj_dict[rxn.reverse_variable] = 2
                        path_opt_vars.append(rxn.forward_variable)
                        path_opt_vars.append(rxn.reverse_variable)

            # 2.2 Create optimization for objective metabolite exchange reaction
            rxn: Reaction = mdl.reactions.get_by_id('SK_'+objective_meta_entry)
            obj_dict[rxn.forward_variable] = 2
            obj_dict[rxn.reverse_variable] = -2

            # 2.3 Create optimization objective for substrate exchange reactions
            # STARTHERE: is this needed?
            substrate_opt_vars = []
            for sub_meta_entry in substrate_meta_entries:
                rxn: Reaction = mdl.reactions.get_by_id('SK_'+sub_meta_entry)
                obj_dict[rxn.forward_variable] = -1
                obj_dict[rxn.reverse_variable] = 1
                substrate_opt_vars.append(rxn.forward_variable)
                substrate_opt_vars.append(rxn.reverse_variable)

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
            for rxn in [x for x in mdl.reactions if x.id != 'SK_'+objective_meta_entry]:
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

            # 2.6 Restricting non-reversible reactions
            # mdl.reactions.get_by_id('R00319').lower_bound = 0
            # mdl.reactions.get_by_id('R03145').lower_bound = 0
            # mdl.reactions.get_by_id('R11074').lower_bound = 0

            # 2.7 Print COBRA model metrics and solve level 1
            print('obj_dict lvl_1',obj_dict)
            mdl.objective = mdl.problem.Objective(Zero, sloppy=True, direction="max")
            mdl.solver.objective.set_linear_coefficients(obj_dict)
            lvl_1_res = self.mass_balance_sln = mdl.optimize()
            print('objective_value lvl_1',lvl_1_res.objective_value)

            fluxes_lvl_1 = self.mass_balance_sln.fluxes
            int_fluxes_lvl_1 = fluxes_lvl_1.loc[[x for x in fluxes_lvl_1.index if 'SK_' not in x]]
            bnd_fluxes_lvl_1 = fluxes_lvl_1.loc[[x for x in fluxes_lvl_1.index if 'SK_' in x]]
            rel_int_fluxes_lvl_1 = int_fluxes_lvl_1.loc[[x for x in int_fluxes_lvl_1.index if round(int_fluxes_lvl_1[x]) != 0]]
            rel_bnd_fluxes_lvl_1 = bnd_fluxes_lvl_1.loc[[x for x in bnd_fluxes_lvl_1.index if round(bnd_fluxes_lvl_1[x]) != 0]]
            print('number of internal rxns used (lvl_1): ',(round(int_fluxes_lvl_1)!=0).sum())
            print('number of boundary rxns used (lvl_1): ',(round(bnd_fluxes_lvl_1)!=0).sum())
            print('int_fluxes (lvl_1): ',rel_int_fluxes_lvl_1.to_markdown())
            print('bnd_fluxes (lvl_1): ',rel_bnd_fluxes_lvl_1.to_markdown())

            # 3. If solution exists, try finding all independent balanced networks
            # 3.1 Restricting to relevant reactions based on previous mass balance solution
            print('starting optimization lvl_2...',rel_int_fluxes_lvl_1.index)
            ### STARTHERE: make second model
            # mdl = Model('mdl')
            rxns_lvl_2 = [x for x in mdl.reactions if x.id in list(rel_int_fluxes_lvl_1.index)]
            # bnd_rxns_lvl_2 = [x for x in mdl.reactions if x.id in list(rel_bnd_fluxes_lvl_1.index)]
            # mdl.add_reactions(rxns_lvl_2+bnd_rxns_lvl_2)
            for rxn in [x for x in mdl.reactions if x not in rxns_lvl_2 and x not in mdl.boundary]:
                mdl.reactions.get_by_id(rxn.id).lower_bound = 0
                mdl.reactions.get_by_id(rxn.id).upper_bound = 0
            print('finished restricting to relevant reactions')

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
                    print('rxn_2',rxn_meta_keys,compound_ids.index(meta.id),self.__cc.get_compound(f"kegg:{meta.id}"),rxn.metabolites[rxn_meta_keys[compound_ids.index(meta.id)]])
                try:
                    eqi_rxn = eq.Reaction(rxn_dict)
                    dGr = abs(self.__cc.dg_prime(eqi_rxn).value.m_as("kJ/mol"))
                    dGr_dict[rxn] = dGr
                    print('dGr',dGr)
                except Exception as e:
                    print('equilibrator warning: ',e)
                    print('abs dGr is assumed negative here')

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
                # obj_dict[new_var] = -10

            # 3.3 preventing the flux of reactions in direction that consume target metabolite
            # target_meta = mdl.metabolites.get_by_id(objective_meta_entry)
            # for rxn in rxns_lvl_2:
            #     stoich = rxn.metabolites
            #     subs = []
            #     prods = []
            #     for key in list(stoich.keys()):
            #         if stoich[key] > 0:
            #             prods.append(key)
            #         elif stoich[key] < 0:
            #             subs.append(key)
            #     if target_meta in prods:
            #         rxn.upper_bound = 0
            #     elif target_meta in subs:
            #         rxn.lower_bound = 0

            print('removing total flux minimization variables and constraints...')
            # mdl.remove_cons_vars([tot_flux_vars,tot_flux_consts])
            for var in tot_flux_vars:
                del obj_dict[var]
            for var in path_opt_vars:
                try:
                    del obj_dict[var]
                except Exception as e:
                    pass
            for var in substrate_opt_vars:
                try:
                    del obj_dict[var]
                except Exception as e:
                    pass

            # 3.4 setting target metabolite flux optimization
            rxn: Reaction = mdl.reactions.get_by_id('SK_'+objective_meta_entry)
            obj_dict[rxn.forward_variable] = len(list(dGr_dict.keys()))
            obj_dict[rxn.reverse_variable] = -len(list(dGr_dict.keys()))

            print('removed total flux minimization variables and constraints')

            # 3.4 Generating optimization target based on enzyme alternative reactions

            # 3.5 Shutting down reactions that are below minimum theoretical stoichiometry
            # max_stoich_val = 0
            # for rxn in rxns_lvl_2:
            #     stoich_vals = list(rxn.metabolites.values())
            #     for val in stoich_vals:
            #         if val > max_stoich_val:
            #             max_stoich_val = val

            # minor_rxn_entries = [x for x in int_fluxes_lvl_1.index if abs(round(int_fluxes_lvl_1[x])) < 1000/max_stoich_val and x in list(map(lambda y: y.id,rxns_lvl_2))]
            # print('minor_rxn_entries',max_stoich_val,minor_rxn_entries)
            # for rxn_entry in minor_rxn_entries:
            #     mdl.reactions.get_by_id(rxn_entry).upper_bound = 0
            #     mdl.reactions.get_by_id(rxn_entry).lower_bound = 0

            # 3.6 Restricts target metabolite production to only the reactions that have the greatest production flux
            direct_target_prod_rxns = {}
            for rxn in rxns_lvl_2:
                stoich = rxn.metabolites
                subs = []
                prods = []
                for key in list(stoich.keys()):
                    if stoich[key] > 0:
                        prods.append(key)
                    elif stoich[key] < 0:
                        subs.append(key)
                if target_meta in prods or target_meta in subs:
                    direct_target_prod_rxns[rxn] = int_fluxes_lvl_1[rxn.id]

            # print('test3',list(direct_target_prod_rxns.values()))
            # max_direct_target_prod_flux = max(list(direct_target_prod_rxns.values()))
            # for rxn in list(direct_target_prod_rxns.keys()):
            #     if direct_target_prod_rxns[rxn] < max_direct_target_prod_flux:
            #         mdl.reactions.get_by_id(rxn.id).upper_bound = 0
            #         mdl.reactions.get_by_id(rxn.id).lower_bound = 0

            # 3.7 Shutting down fluxes that produce substrate metabolites
            # substrate_metas = [mdl.metabolites.get_by_id(x) for x in substrate_meta_entries]
            # for rxn in rxns_lvl_2:
            #     stoich = rxn.metabolites
            #     subs = []
            #     prods = []
            #     for key in list(stoich.keys()):
            #         if stoich[key] > 0:
            #             prods.append(key)
            #         elif stoich[key] < 0:
            #             subs.append(key)
            #     for meta in substrate_metas:
            #         print('subs and prods',meta,subs,prods)
            #         if meta in prods:
            #             mdl.reactions.get_by_id(rxn.id).upper_bound = 0
            #         elif meta in subs:
            #             mdl.reactions.get_by_id(rxn.id).lower_bound = 0

            # 3.8 Print COBRA model metrics and solve level 2
            print('obj_dict lvl_2',obj_dict)
            mdl.objective = mdl.problem.Objective(Zero, sloppy=True, direction="max")
            mdl.solver.objective.set_linear_coefficients(obj_dict)
            lvl_2_res = self.mass_balance_sln = mdl.optimize()
            print('objective_value lvl_2',lvl_2_res.objective_value)

            fluxes_lvl_2 = self.mass_balance_sln.fluxes
            int_fluxes_lvl_2 = fluxes_lvl_2.loc[[x for x in fluxes_lvl_2.index if 'SK_' not in x]]
            bnd_fluxes_lvl_2 = fluxes_lvl_2.loc[[x for x in fluxes_lvl_2.index if 'SK_' in x]]
            rel_int_fluxes_lvl_2 = int_fluxes_lvl_2.loc[[x for x in int_fluxes_lvl_2.index if round(int_fluxes_lvl_2[x]) != 0]]
            rel_bnd_fluxes_lvl_2 = bnd_fluxes_lvl_2.loc[[x for x in bnd_fluxes_lvl_2.index if round(bnd_fluxes_lvl_2[x]) != 0]]
            print('number of internal rxns used (lvl_2): ',(round(int_fluxes_lvl_2)!=0).sum())
            print('number of boundary rxns used (lvl_2): ',(round(bnd_fluxes_lvl_2)!=0).sum())
            print('int_fluxes (lvl_2): ',rel_int_fluxes_lvl_2.to_markdown())
            print('bnd_fluxes (lvl_2): ',rel_bnd_fluxes_lvl_2.to_markdown())

            print('starting optimization lvl_3...',rel_int_fluxes_lvl_2.index)
            rxns_lvl_3 = [x for x in mdl.reactions if x.id in list(rel_int_fluxes_lvl_2.index)]
            for rxn in [x for x in mdl.reactions if x not in rxns_lvl_3 and x not in mdl.boundary]:
                mdl.reactions.get_by_id(rxn.id).lower_bound = 0
                mdl.reactions.get_by_id(rxn.id).upper_bound = 0
            print('finished restricting to relevant reactions')

            # 4. If small gas metas are being produced and consumed
            # small_gas_metas_circ = {}
            # small_gas_metas = [mdl.metabolites.get_by_id(x) for x in self.__small_gas_metas]
            # # 4.1 Prevents production of small gas metabolites if exchange is producing, vice versa
            # small_gas_meta_exchange_dir = {}
            # print('test3',rel_bnd_fluxes_lvl_2.index)
            # print('test2',rel_bnd_fluxes_lvl_2.index.values)
            # print('test',[x for x in small_gas_metas if 'SK_'+x.id in list(rel_bnd_fluxes_lvl_2.index.values)])
            # for small_meta in [x for x in small_gas_metas if 'SK_'+x.id in list(rel_bnd_fluxes_lvl_2.index.values)]:
            #     if rel_bnd_fluxes_lvl_2.loc['SK_'+small_meta.id] > 0:
            #         small_gas_meta_exchange_dir[small_meta.id] = 'out'
            #     elif rel_bnd_fluxes_lvl_2.loc['SK_'+small_meta.id] < 0:
            #         small_gas_meta_exchange_dir[small_meta.id] = 'in'
            # print('rxns that are edited',small_gas_meta_exchange_dir,[x for x in rxns_lvl_2 if 'SK_' not in x.id])
            # for rxn in [x for x in rxns_lvl_2 if 'SK_' not in x.id]:
            #     stoich = rxn.metabolites
            #     subs = []
            #     prods = []
            #     print('stoich',list(stoich.keys()))
            #     for key in list(stoich.keys()):
            #         if stoich[key] > 0:
            #             prods.append(key.id)
            #         elif stoich[key] < 0:
            #             subs.append(key.id)
            #     for meta in [x for x in small_gas_metas if x.id in list(small_gas_meta_exchange_dir.keys())]:
            #         print('rxn.id',meta.id,rxn.id,subs,prods,small_gas_meta_exchange_dir[meta.id])
            #         if rxn.id not in list(small_gas_metas_circ.keys()):
            #             small_gas_metas_circ[rxn.id] = []
            #         if small_gas_meta_exchange_dir[meta.id] == 'in':
            #             if meta.id in prods:
            #                 small_gas_metas_circ[rxn.id].append('upper')
            #             elif meta.id in subs:
            #                 small_gas_metas_circ[rxn.id].append('lower')
            #         elif small_gas_meta_exchange_dir[meta.id] == 'out':
            #             if meta.id in prods:
            #                 small_gas_metas_circ[rxn.id].append('lower')
            #             elif meta.id in subs:
            #                 small_gas_metas_circ[rxn.id].append('upper')

            # print('restricting fluxes for small gas metas',small_gas_metas_circ)
            # if len(list(small_gas_metas_circ.keys())) > 0:
            #     for key in list(small_gas_metas_circ.keys()):
            #         if 'upper' in small_gas_metas_circ[key]:
            #             mdl.reactions.get_by_id(key).lower_bound = 0
            #         elif 'lower' in small_gas_metas_circ[key]:
            #             mdl.reactions.get_by_id(key).upper_bound = 0

            # 4.1 choose the production pathway of target meta that generates the highest yield
            max_target_meta_prod_rxn = {}
            for rxn in [x for x in rxns_lvl_3 if 'SK_' not in x.id]:
                if objective_meta_entry in list(map(lambda y: y.id,list(rxn.metabolites.keys()))):
                    max_target_meta_prod_rxn[rxn.id] = abs(rel_int_fluxes_lvl_2.loc[rxn.id])
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

                # obj_dict[mdl.reactions.get_by_id(max_flux_rxn_id).forward_variable] = len(list(dGr_dict.keys()))
                # obj_dict[mdl.reactions.get_by_id(max_flux_rxn_id).reverse_variable] = -len(list(dGr_dict.keys()))

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

        # Compare balanced networks that get to target metabolite from each given substrate, try to consolidate into most efficient 
        # pathway that uses some combination of all of the substrate metabolites

        return

    # balances reaction network
    def __rank_candidate_networks(self,objective_meta_entries:list[str],
                                         substrate_meta_entries:list[str],
                                         boundary_metabolite_entries:list[str],
                                         manual_zero_flux_rxns:list[str],
                                         opt_weights:dict[str,float]) -> None:
        """STARTHERE: try methods to get better network:
        1. Prune (set to 0) the smallest fluxes
        2. Brute-force required substrate combinations
        3. Try assuming no other exchange other than those outputted
        4. Try assuming that all substrates are utilized
        """

        # 1. Set constraints for COBRA model - fully shut down certain reaction fluxes


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

        # find potential paths between substrate(s) and objective metabolite
        candidate_paths = self.__find_shortest_simple_paths(objective_meta_entry=objective_metabolite_entry,
                                                            substrate_meta_entries=substrate_metabolite_entries,
                                                            min_enzyme_ct=min_enzyme_ct,
                                                            max_enzyme_ct=max_enzyme_ct)

        # balance all candidate paths
        self.__balance_candidate_networks(candidate_paths=candidate_paths,
                                          objective_meta_entry=objective_metabolite_entry,
                                          substrate_meta_entries=substrate_metabolite_entries)

        # self.__rank_candidate_networks(objective_meta_entries=objective_metabolite_entries,
        #                                       substrate_meta_entries=substrate_metabolite_entries,
        #                                       boundary_metabolite_entries=objective_metabolite_entries+substrate_metabolite_entries,
        #                                       manual_zero_flux_rxns=[],
        #                                       opt_weights=optimization_weights)

        return self.__COBRA_model

    def assess_enzyme_promiscuity(self) -> None:
        return