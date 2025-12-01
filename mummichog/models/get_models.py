'''
metabolic models from
a. Python import here (metabolicModels.py)
b. JSON files here under json/
c. external database JMS/Azimuth DB (future)

from jms.modelConvert import convert_json_model
from jms.empiricalCpds import load_epds_from_json

>>> from _metabolic_models import *
>>> models.keys()
dict_keys(['Staphylococcus_epidermidis_ATCC_12228_massInferred', 'metabolicModel_EBI_20210602_Neisseria_meningitidis_alpha14', 'metabolicModel_AGORA_20210512_Escherichia_coli_O157_H7_str_Sakai', 'metabolicModel_EBI_20210602_Propionibacterium_acnes_SK137', 'metabolicModel_az_HumanGEM_20220302_noCompartmentalization', 'worm_model_icel1273'])

>>> for m in models:
...   print(models[m].keys())
... 
dict_keys(['id', 'version', 'Compounds', 'dict_cpds_def', 'metabolic_rxns', 'cpd_edges', 'edge2rxn', 'edge2enzyme', 'metabolic_pathways', 'cpd2pathways'])
dict_keys(['id', 'version', 'Compounds', 'dict_cpds_def', 'metabolic_rxns', 'cpd_edges', 'edge2rxn', 'edge2enzyme', 'metabolic_pathways', 'cpd2pathways'])
dict_keys(['id', 'version', 'Compounds', 'dict_cpds_def', 'metabolic_rxns', 'cpd_edges', 'edge2rxn', 'edge2enzyme', 'metabolic_pathways', 'cpd2pathways'])
dict_keys(['id', 'version', 'Compounds', 'dict_cpds_def', 'metabolic_rxns', 'cpd_edges', 'edge2rxn', 'edge2enzyme', 'metabolic_pathways', 'cpd2pathways'])
dict_keys(['id', 'version', 'Compounds', 'dict_cpds_def', 'metabolic_rxns', 'cpd_edges', 'edge2rxn', 'edge2enzyme', 'metabolic_pathways', 'cpd2pathways'])
dict_keys(['metabolic_rxns', 'cpd_edges', 'metabolic_pathways', 'Compounds', 'dict_cpds_def', 'cpd2pathways', 'edge2enzyme', 'edge2rxn', 'version', 'dict_cpds_mass'])
>>> 
>>> 
>>> 
>>> list(models['metabolicModel_EBI_20210602_Neisseria_meningitidis_alpha14']['cpd2pathways'].items())[:10]
[]
>>> list(models['metabolicModel_EBI_20210602_Neisseria_meningitidis_alpha14']['edge2enzyme'].items())[:10]
[]
>>> models['worm_model_icel1273']['metabolic_pathways'][:2]
[{'cpds': ['C00267', 'C00221', 'C00031'], 'rxns': ['RCC0138', 'RCC0137'], 'ecs': ['NA'], 'name': '', 'id': 'iCEL1273pathway1'}, {'cpds': ['C00025', 'C00024', 'C00148', 'C00026', 'M00016', 'C00007', 'C00006', 'C00005', 'C00004', 'C00003', 'C00027', 'C00001', 'C03912', 'C02946', 'M00148', 'C00033', 'M01165', 'C00437', 'C05947', 'C05946', 'M01352', 'M03287', 'C00334', 'C00014', 'C00010', 'M00077', 'C02714', 'M00001', 'C00011', 'M00003', 'M00002', 'M00005', 'M00004', 'M00006', 'M00009', 'M00008', 'C00134', 'C00077', 'M00025', 'M00026', 'M03912', 'M05947', 'M05946', 'M05938', 'M04281', 'M01157', 'C04281', 'M00080', 'C00080', 'C00042', 'C01157', 'C05936'], 'rxns': ['RC05050', 'RC01252', 'RC05052', 'RC01251', 'RC03293', 'RC01987', 'RMC0016', 'RC03291', 'RMC0015', 'RM05052', 'RM00239', 'RC00670', 'RM01248', 'RM00667', 'RC01248', 'RC04025', 'RC01154', 'RC00669', 'RM04444', 'RM04445', 'RM03293', 'RM03291', 'RM00245', 'RM05051', 'RM03313', 'RM01251', 'RM01253'], 'ecs': ['1.2.1.88', '1.5.99.8', '2.6.1.13', '1.14.11.2', '3.5.1.14', 'spontaneous', '1.5.1.2', '1.4.3.4', '1.2.1.41', '2.6.1.1', '4.1.1.17', '2.7.2.11', '1.5.1.12', '3.5.1.63', '1.5.-.-', '1.2.1.3', '3.5.1.16', '2.3.1.57'], 'name': 'Arginine and proline metabolism', 'id': 'iCEL1273pathway7'}]
>>> 
>>> models['metabolicModel_EBI_20210602_Neisseria_meningitidis_alpha14']['metabolic_pathways'][:2]
[]
>>> models['metabolicModel_EBI_20210602_Neisseria_meningitidis_alpha14']['metabolic_rxns'][:2]
[{'id': 'MNXR2184_i', 'reactants': ['bigg_fdp_i'], 'products': ['MNXM74_i', 'bigg_dhap_i']}, {'id': 'MNXR648_i', 'reactants': ['bigg_h_i', 'bigg_nadh_i', 'bigg_acald_i'], 'products': ['bigg_etoh_i', 'bigg_nad_i']}]


'''

# import json
import networkx as nx
from .metabolicModels import metabolicModels

def get_remote_metabolic_model(db=None, model_id=None):
    '''
    pull model from the remote db
    '''
    pass

def get_metabolic_model(model='human_model_mfn'):
    '''
    To-do:
    handling models from JMS and other sources
    
    '''
    return metabolicNetwork(metabolicModels[ model ])


class metabolicNetwork:
    '''
    Metabolite-centric metabolic model 
    Theoretical model, not containing user data

    This is from # from metDataModel.mummichog import metabolicNetwork


    '''
    def __init__(self, MetabolicModel):
        '''
        Initiation of metabolic network model.
        Building Compound index.
        Parsing input files.
        Matching m/z - Compound.
        
        MetabolicModel['Compounds'] are subset of cpds in network/pathways with mw.
        Not all in total_cpd_list has mw.
        '''
        #print_and_loginfo( "Loading metabolic network %s..." %MetabolicModel.version ) # version from metabolic model
        
        self.MetabolicModel = MetabolicModel
        self.network = self.build_network(MetabolicModel['cpd_edges'])
        
        self.version = MetabolicModel['version']
        self.Compounds = MetabolicModel['Compounds']
        self.metabolic_pathways = MetabolicModel['metabolic_pathways']
        self.dict_cpds_def = MetabolicModel['dict_cpds_def']
        self.cpd2pathways = MetabolicModel['cpd2pathways']
        self.edge2enzyme = MetabolicModel['edge2enzyme']
        self.total_cpd_list = self.network.nodes()
        
        
    def build_network(self, edges):
        return nx.from_edgelist( edges )
        

    def get_pathways(self):
        pass

