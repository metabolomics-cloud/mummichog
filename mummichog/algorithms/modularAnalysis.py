#
# --------------------------------------------------------
#
# module analysis
#
import sys
import logging
import random
import itertools
import numpy as np
from scipy import stats
import networkx as nx

SEARCH_STEPS = 4
MODULE_SIZE_LIMIT = 100
SIGNIFICANCE_CUTOFF = 0.05   # to get from parameters later
USE_DEBUG = False

def find_communities(G, method='louvain'):
    '''
    Wrapper of nx.community functions
    
    '''
    if method == 'louvain':
        return nx.community.louvain_communities(G)
    elif method == 'Clauset-Newman-Moore':
        return nx.community.greedy_modularity_communities(G)
    elif method == 'girvan_newman':
        # no optiaml, returning all levels
        coms = nx.community.girvan_newman(G)
        return [list(x) for x in coms]
    else:
        print("Warning, method not implemented")
        


class Mmodule:
    '''
    Metabolites by their connection in metabolic network.
    A module is a subgraph, while modularity is calculated in 
    the background of reference hsanet.
    
    
    need to record sig EmpCpds
    
    '''
    def __init__(self, network, subgraph, TrioList):
        '''
        TrioList (seeds) format: [(M.row_number, EmpiricalCompounds, Cpd), ...]
        to keep tracking of where the EmpCpd came from (mzFeature).
        
        network is the total parent metabolic network
        '''
        self.network = network
        self.num_ref_edges = self.network.number_of_edges()
        self.num_ref_nodes = self.network.number_of_nodes()
        self.graph = subgraph.copy()
        
        seed_cpds = [x[2] for x in TrioList]
        self.shave(seed_cpds)
        self.nodestr = self.make_nodestr()
        self.N_seeds = len(seed_cpds)
        self.A = self.activity_score(seed_cpds, self.get_num_EmpCpd(TrioList))
    
    def activity_score(self, seed_cpds, num_EmpCpd):
        '''
        A * (Ns/Nm)
        A = Newman-Girvan modularity score
        Ns = number of input cpds in module M
        Nm = number of total cpds in M
        Ns/Nm can be corrected as (Ns/total input size)/(Nm/network size), however,
        this normalization factor holds the same in permutations. 
        Use 100 here for network size/total input size.
        
        To reduce bias towards larger modules in Q:
        np.sqrt(len(seed_cpds)/Nm) * 
        
        Ns is now controlled by number of empiricalCompounds
        '''
        Ns = num_EmpCpd
        Nm = float(self.graph.number_of_nodes())
        if Nm > 0:
            self.compute_modularity()
            return np.sqrt(self.N_seeds/Nm) *self.Q * (Ns/Nm) * 100
        else:
            return 0
        
        
    def get_num_EmpCpd(self, TrioList):
        new = []
        subgraph_nodes = self.graph.nodes()
        for x in TrioList:
            if x[2] in subgraph_nodes:
                new.append(x[1])
                
        return len(set(new))
        
        
    def compute_modularity(self):
        '''
        To compute Newman-Girvan modularity for a single module,
        in reference to the whole network.
        '''
        m = self.num_ref_edges
        Nodes = self.graph.nodes()
        expected = 0
        for ii in Nodes:
            for jj in Nodes:
                if ii != jj:
                    expected += self.network.degree(ii) * self.network.degree(jj)
                    
        expected /= (4.0 * m)
        self.Q = (self.graph.number_of_edges() - expected) / m
    
    def test_compute_modularity(self):
        '''
        Alternative modularity measure as 
        edges in module over all edges on the same nodes
        '''
        m = float(self.graph.number_of_edges())
        expected = 0
        for ii in self.graph.nodes(): expected += self.network.degree(ii)
        self.Q = 2 * m * (np.sqrt(self.graph.number_of_nodes())) / expected
        
        
    def shave(self, seed_cpds):
        '''
        shave off nodes that do not connect seeds, i.e.
        any node with degree = 1 and is not a seed, iteratively.
        '''
        nonseeds = [x for x in self.graph.nodes() if x not in seed_cpds]
        excessive = [x for x in nonseeds if self.graph.degree(x)==1]
        while excessive:
            for x in excessive: 
                self.graph.remove_node(x)
            nonseeds = [x for x in self.graph.nodes() if x not in seed_cpds]
            excessive = [x for x in nonseeds if self.graph.degree(x)==1]


    def make_nodestr(self):
        '''
        create an identifier using nodes in sorted order
        '''
        Nodes = list(self.graph.nodes)
        Nodes.sort()
        return ''.join(Nodes)

    def export_network_txt(self, met_model, filename):
        '''
        To use .txt for Cytoscape 3, no need for .sif any more.
        Edges are strings now as switching to JSON compatible.
        '''
        s = 'SOURCE\tTARGET\tENZYMES\n'
        for e in self.graph.edges():
            s += e[0] + '\t' + e[1] + '\t' + met_model.edge2enzyme.get(','.join(sorted(e)), '') + '\n'
        
        out = open(filename, 'w')
        out.write(s)
        out.close()


        
class ModularAnalysis:
    '''
    1) Find modules from input list by connecting paths < 4;
    compute activity score that combines modularity and enrichment.
    2) Permutations by randomly selecting features from ref_mzlist;
    compute p-values based on permutations.
    
    Module analysis will still be in the compound space, as network model is defined by compound edges.
    Need better tracking the mapping btw compound and EmpiricalCompounds
    
    Tested to add a generator from EmpiricalCompounds to a bunch of combinations, 
    i.e. only one cpd from Ecpd is used at a time towards module analysis
    But it's too slow to be practical.
    
    ...

    Attributes
    ----------

    Methods
    -------
    
    '''
    def __init__(self, mixedNetwork):
        '''
        mapping btw (mzfeature, cpd) has to be via ListOfEmpiricalCompounds, 
        so that cpd can be tracked back to EmpiricalCompounds
        
        
        '''
        self.mixedNetwork = mixedNetwork
        self.network = mixedNetwork.model.network
        self.paradict = mixedNetwork.data.paradict
        
        # both using row_numbers
        self.ref_featurelist = self.mixedNetwork.features
        self.significant_features = self.mixedNetwork.significant_features
        self.significant_Trios = self.mixedNetwork.TrioList
        

    def dispatch(self):
        '''
        Only real modules are saved in total. 
        Permutated modules are not saved but their scores are recorded. 
        
        
        print("what found? -", self.modules_from_significant_features[3].nodestr)
        
        
        '''
        s = "\nModular Analysis, using %d permutations ..." %self.paradict['permutation']
        print (s)
        self.modules_from_significant_features = self.run_analysis_real()
        self.permuation_mscores = self.do_permutations(self.paradict['permutation'])

        self.rank_significance()        
        #for M in self.top_modules: print(M, M.A, nx.average_node_connectivity(M.graph))


    def run_analysis_real(self):
        return self.find_modules( self.significant_Trios )
        

    def do_permutations(self, num_perm):
        '''
        Run num_perm permutations on ref featurelist;
        populate activity scores from random modules in self.permuation_mscores
        
        '''
        permuation_mscores = []
        N = len(self.significant_features)
        for ii in range(num_perm):
            sys.stdout.write( ' ' + str(ii+1))
            sys.stdout.flush()
            random_trios = self.mixedNetwork.batch_rowindex_EmpCpd_Cpd( 
                                            random.sample(self.ref_featurelist, N) )
            permuation_mscores += [x.A for x in self.find_modules(random_trios)] or [0]
            
        return permuation_mscores
            
    def __generator_EmpiricalCompounds_cpds__(self, Ecpds):
        '''
        return the combinations of one cpd drawn from each EmpiricalCompound, N = len(Ecpds)
        as iterator
        This function is not used now, because too many combinations make software too slow.
        '''
        return itertools.product(*[ E.compounds for E in Ecpds ])


    def find_modules(self, TrioList):
        '''
        get connected nodes in up to 4 steps.
        modules are set of connected subgraphs plus split moduels within.
        A shaving step is applied to remove excessive nodes that do not 
        connect seeds (thus Mmodule initiation may reduce graph size). 
        A module is only counted if it contains more than one seeds.
        
        TrioList format: [(M.row_number, EmpiricalCompounds, Cpd), ...]
        '''
        global SEARCH_STEPS, MODULE_SIZE_LIMIT
        seeds = [x[2] for x in TrioList]      # use cpd space
        modules, modules2, module_nodes_list = [], [], []

        for ii in range(SEARCH_STEPS):
            edges = nx.edges(self.network, seeds)
            if ii == 0:
                # step 0, counting edges connecting seeds
                edges = [x for x in edges if x[0] in seeds and x[1] in seeds]
                new_network = nx.from_edgelist(edges)
                
            else:
                # step 1, 2, 3, ... growing to include extra steps/connections
                new_network = nx.from_edgelist(edges)
                seeds = new_network.nodes()
            
            for sub in nx.connected_components(new_network):
                sub = new_network.subgraph(sub)
                if 3 < sub.number_of_nodes() < MODULE_SIZE_LIMIT:
                    M = Mmodule(self.network, sub, TrioList)
                    modules.append(M)
                
        # add modules split from modules
        if USE_DEBUG:
            logging.info( '# initialized network size = %d' %len(seeds) )
            # need export modules for comparison to heinz
            self.__export_debug_modules__( modules )
            
        for sub in modules:
            if sub.graph.number_of_nodes() > 5:
                modules2 += [Mmodule(self.network, x, TrioList)
                             for x in self.__split_modules__(sub.graph)]
        
        new = []
        for M in modules + modules2:
            if M.graph.number_of_nodes() > 3 and M.nodestr not in module_nodes_list:
                new.append(M)
                module_nodes_list.append(M.nodestr)
                if USE_DEBUG: logging.info( str(M.graph.number_of_nodes()) + ', ' + str(M.A) )

        return new


    def __export_debug_modules__(self, modules):
        '''
        write out initial modules, to be split by alternative algorithm
        '''
        s = ''
        for M in modules: s += M.make_sif_str()
        out = open( os.path.join(self.modules_dir, 'debug_modules.sif'), 'a' )
        out.write(s + '#\n')
        out.close()
        
    def __split_modules__(self, g):
        '''
        return nx.graph instance after splitting the input graph
        by Newman's spectral split method
        Only modules more than 3 nodes are considered as good small modules 
        should have been generated in 1st connecting step.
        '''
        return [nx.subgraph(g, x) for x in find_communities(g) if len(x) > 3]


    def rank_significance(self):
        '''
        compute p-values of modules. Either model based:
        scores of random modules are fitted to a Gamma distribution,
        p-value is calculated from CDF.
        Or rank based.
        '''
        print("\nNull distribution is estimated on %d random modules" 
                          %len(self.permuation_mscores))
        print("User data yield %d network modules" 
                          %len(self.modules_from_significant_features))
        
        if self.paradict['modeling'] == 'gamma':
            a, loc, scale = stats.gamma.fit(self.permuation_mscores)
            if USE_DEBUG:
                logging.info( 'Gamma fit parameters a, loc, scale = ' + ', '.join([str(x) for x in [a, loc, scale]]) )
            
            for M in self.modules_from_significant_features:
                M.p_value = 1 - stats.gamma.cdf(M.A, a, loc, scale)
            
        else:
            for M in self.modules_from_significant_features:
                M.p_value = self.__calculate_p__(M.A, self.permuation_mscores)
                
        top_modules = [M for M in self.modules_from_significant_features if M.p_value < SIGNIFICANCE_CUTOFF]
        self.top_modules = sorted(top_modules, key=lambda M: M.p_value)
        
        
        


    def __calculate_p__(self, x, record):
        '''
        calculate p-value based on the rank in record of scores
        '''
        total_scores = [x] + record
        total_scores.sort(reverse=True)
        D = len(record) + 1.0
        return (total_scores.index(x)+1)/D
  
    def collect_hit_Trios(self):
        '''
        get [(mzFeature, EmpiricalCompound, cpd),...] for top_modules.
        Update EmpCpd chosen compounds.
        '''
        overlap_Cpds = []
        for M in self.top_modules:
            overlap_Cpds += M.graph.nodes()
        
        overlap_Cpds = set(overlap_Cpds)
        new = []
        for T in self.significant_Trios:
            if T[2] in overlap_Cpds:
                # yet to sort this out
                
                # T[1].update_chosen_cpds(T[2])
                # T[1].designate_face_cpd()
                new.append(T)
        
        return new



    def to_json(self):
        '''
        JSON export to be consumed by downstream functions


        '''
        L = []
        for M in self.top_modules:
            L.append({
                'p-value': M.p_value,
                'edges': [','.join(e) for e in M.graph.edges()],
                'nodes': list(M.graph.nodes)
            })
        return L
