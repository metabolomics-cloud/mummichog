#
# --------------------------------------------------------
#
# activity network analysis
#

import networkx as nx


class ActivityNetwork:
    '''
    Tally cpds responsible for significant pathways and modules,
    and build an acitvity network to represent top story in the data.
    Remove singletons in network. Try no more than 3 steps. AN can get too big or too small.
    
    '''
    def __init__(self, mixedNetwork, hit_Trios):
        '''
        Build a consensus network for hit_Trios,
        [(mzFeature, EmpiricalCompound, cpd),...] for top_modules and sig pathways.
        hit_Trios = set(PA.collect_hit_Trios() + MA.collect_hit_Trios())
        '''
        # also update to mixedNetwork
        mixedNetwork.hit_Trios = hit_Trios
        self.mixedNetwork = mixedNetwork
        self.network = mixedNetwork.model.network
        
        nodes = [x[2] for x in hit_Trios]
        self.activity_network = self.build_activity_network(nodes)

    def build_activity_network(self, nodes, cutoff_ave_conn = 0.5, expected_size = 10):
        '''
        Get a network with good connections in no more than 3 steps.
        No modularity requirement for 1 step connected nodes.
        '''
        an = nx.subgraph(self.mixedNetwork.model.network, nodes)
        if nodes:
            sub1 = self.__get_largest_subgraph__(an)
            if sub1.number_of_nodes() > expected_size:
                print("\nActivity network was connected in 1 step.")
                return sub1
            
            else:   # expand 1 or 2 steps
                edges = nx.edges(self.mixedNetwork.model.network, nodes)
                new_network = self.__get_largest_subgraph__( nx.from_edgelist(edges) )
                conn = self.__get_ave_connections__(new_network)
                if an.number_of_nodes() > MODULE_SIZE_LIMIT or conn > cutoff_ave_conn:
                    print("\nActivity network was connected in 2 steps.")
                    return new_network
                else:
                    edges = nx.edges(self.mixedNetwork.model.network, new_network.nodes())
                    new_network = self.__get_largest_subgraph__( nx.from_edgelist(edges) )
                    conn = self.__get_ave_connections__(new_network)
                    if conn > cutoff_ave_conn:
                        print("\nActivity network was connected in 3 steps.")
                        return new_network
                    else:
                        return an
        else:
            return an
        
    def export_network_txt(self, met_model, filename):
        s = 'SOURCE\tTARGET\tENZYMES\n'
        for e in self.activity_network.edges():
            s += e[0] + '\t' + e[1] + '\t' + met_model.edge2enzyme.get(','.join(sorted(e)), '') + '\n'
        
        out = open(filename, 'w')
        out.write(s)
        out.close()
        
    def __get_largest_subgraph__(self, an):
        '''
        connected_component_subgraphs likely to return sorted subgraphs. Just to be sure here.
        '''
        connected_component_subgraphs = (an.subgraph(c) for c in nx.connected_components(an))
        return max(connected_component_subgraphs, key=len)
        
    def __get_ave_connections__(self, N):
        '''
        nx.average_node_connectivity(G) is too slow; use self.__get_ave_connections__()
        '''
        try:                #Avoid ZeroDivisionError
            return N.number_of_edges()/float(N.number_of_nodes())
        except:
            return 0


    
    def to_json(self):
        '''
        Convert result to dataframes, easy JSON export to be consumed by downstream functions


        '''

        return  [','.join(e) for e in self.activity_network.edges()]
        
