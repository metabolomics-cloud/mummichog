'''
'''
import sys
from scipy import stats


class metabolicPathway:
    def __init__(self):
        self.id = ''
        self.name = ''
        self.rxns = []
        self.ecs = []
        self.ec_num = 0
        self.cpds = []
        self.cpd_num = 0
        
        self.selected_features = []
        self.combined_pvalue = 0    # -log10()
        
    def str_import(self, s):
        '''
        function to import pathway from plain text;
        excluding currency metabolites.
        Not used for now.
        '''
        a = s.rstrip().split('\t')
        self.id = a[0]
        self.name = a[1]
        self.rxns = a[2].split(';')
        self.ecs = a[3].split(';')
        self.ec_num = int(a[4])
        cpds = a[5].split(';')
        self.cpds = [x for x in cpds if x not in currency]
        self.cpd_num = len(self.cpds)
        
    def json_import(self, j):
        '''
        function to import pathway from JSON format.
        '''
        self.id = j['id']
        self.name = j['name']
        self.rxns = j['rxns']
        self.ecs = j['ecs']
        self.ec_num = len(j['ecs'])
        self.cpds = j['cpds']
        self.cpd_num = len(j['cpds'])


class PathwayAnalysis:
    '''
    Pathway enrichment analysis, considering uncertainty in metabolite annotation.
    
    p-value is from Fisher exact test, 
    adjusted by resampling method in 
    GF Berriz, OD King, B Bryant, C Sander & FP Roth. 
    Characterizing gene sets with FuncAssociate. 
    Bioinformatics 19(18):2502-2504 (2003)
    
    '''
    def __init__(self, pathways, mixedNetwork):
        '''
        mixedNetwork contains both user input data, metabolic model,
        and mapping btw (mzFeature, EmpiricalCompound, cpd)
        
        '''
        self.mixedNetwork = mixedNetwork
        self.network = mixedNetwork.model.network
        self.paradict = mixedNetwork.data.paradict
        
        self.pathways = self.get_pathways(pathways)
        self.resultListOfPathways = []          # will store result of pathway analysis
        
        # to help track wehre sig cpd comes from
        self.TrioList = self.mixedNetwork.TrioList
        self.significant_EmpiricalCompounds = set([x[1] for x in self.TrioList])
        
        self.ListOfEmpiricalCompounds = mixedNetwork.ListOfEmpiricalCompounds
        self.total_number_EmpiricalCompounds = len(self.ListOfEmpiricalCompounds)

        print("\nPathway Analysis...")
        
        
    def get_pathways(self, pathways):
        '''
        convert pathways in JSON formats (import from .py) to list of Pathway class.
        Adding list of EmpiricalCompounds per pathway, which reflects the measured pathway coverage.
        '''
        new = []
        for j in pathways:
            P = metabolicPathway()
            P.json_import(j)
            P.adjusted_p = ''
            P.EmpiricalCompounds = self.__get_empiricalCompounds_by_cpds__(P.cpds)
            new.append(P)
        return new
        

    def __get_empiricalCompounds_by_cpds__(self, cpds):
        '''
        Mapping cpds to empirical_cpds. Also used for counting EmpCpds for each Pathway.
        '''
        cpds_empirical = []
        for c in cpds: cpds_empirical += self.mixedNetwork.Compound_to_EmpiricalCompounds.get(c, [])
        return set(cpds_empirical)
        
        
    def do_permutations(self, pathways, num_perm):
        '''
        Modified from Berriz et al 2003 method.
        After collecting p-values from resampling, do a Gamma fit.
        
        Permutation is simplified in version 2; no more new TableFeatures instances.
        
        
        May consider fitting Gamma at log scale, to be more accurate --
        
        '''
        self.permutation_record = []
        print("Resampling, %d permutations to estimate background ..." 
                          %num_perm)
        
        # this is feature number not cpd number
        N = len(self.mixedNetwork.significant_features)
        for ii in range(num_perm):
            sys.stdout.write( ' ' + str(ii + 1))
            sys.stdout.flush()
            random_Trios = self.mixedNetwork.batch_rowindex_EmpCpd_Cpd( 
                                random.sample(self.mixedNetwork.mzrows, N) )
            query_EmpiricalCompounds = set([x[1] for x in random_Trios])
            self.permutation_record += (self.__calculate_p_ermutation_value__(
                query_EmpiricalCompounds, pathways))
        
        print("\nPathway background is estimated on %d random pathway values" 
                          %len(self.permutation_record))
        


    def __calculate_p_ermutation_value__(self, query_EmpiricalCompounds, pathways):
        '''
        calculate the FET p-value for all pathways.
        But not save anything to Pathway instances.
        '''
        p_of_pathways = [ ]
        query_set_size = len(query_EmpiricalCompounds)
        total_feature_num = self.total_number_EmpiricalCompounds
        
        for P in pathways:
            overlap_features = query_EmpiricalCompounds.intersection(P.EmpiricalCompounds)
            overlap_size = len(overlap_features)
            ecpd_num = len(P.EmpiricalCompounds)
            if overlap_size > 0:
                negneg = total_feature_num + overlap_size - ecpd_num - query_set_size
                p_val = stats.fisher_exact([[overlap_size, query_set_size - overlap_size],
                                       [ecpd_num - overlap_size, negneg]], 'greater')[1]
                p_of_pathways.append(p_val)
            else: 
                p_of_pathways.append(1)
                
        return p_of_pathways


    def get_adjust_p_by_permutations(self, pathways):
        '''
        EASE score is used as a basis for adjusted p-values,
        as mummichog encourages bias towards more hits/pathway.
        pathways were already updated by first round of Fisher exact test,
        to avoid redundant calculations.
        "Adjusted_p" is not an accurate term. It's rather a permutation based empirical p-value.
        '''
        self.do_permutations(pathways, self.paradict['permutation'])
        
        if self.paradict['modeling'] == 'gamma':
            #vector_to_fit = [-np.log10(x) for x in self.permutation_record if x < 1]
            vector_to_fit = -np.log10(np.array(self.permutation_record))
            self.gamma = stats.gamma.fit(vector_to_fit)
            a, loc, scale = self.gamma
            
            for P in pathways: 
                P.adjusted_p = self.__calculate_gamma_p__(a, loc, scale, P.p_EASE)
        else:
            for P in pathways: P.adjusted_p = self.__calculate_p__(P.p_EASE, self.permutation_record)
        return pathways
        

    def __calculate_p__(self, x, record):
        '''
        calculate p-value based on the rank in record of permutation p-values
        '''
        total_scores = [x] + record
        total_scores.sort()
        D = len(record) + 1.0
        return (total_scores.index(x)+1)/D
    
    def __calculate_gamma_p__(self, a, loc, scale, x):
        '''
        Use -log10 scale for model fitting
        '''
        return 1 - stats.gamma.cdf(-np.log10(x), a, loc, scale)
    
    
    def cpd_enrich_test(self):
        '''
        Fisher Exact Test in cpd space, after correction of detected cpds.
        Fisher exact test is using scipy.stats.fisher_exact
        for right-tail p-value:
        >>> stats.fisher_exact([[12, 5], [29, 2]], 'greater')[1]
        0.99452520602188932
        
        query size is now counted by EmpiricalCompounds.
        adjusted_p should be model p-value, not fdr.
        This returns a list of Pathway instances, with p-values.
        
                        P.p_EASE = stats.fisher_exact([[max(0, overlap_size - 1), query_set_size - overlap_size],
                                   [ecpd_num - overlap_size + 1, negneg]], 'greater')[1]
        '''
        FET_tested_pathways = []
        qset = self.significant_EmpiricalCompounds
        query_set_size = len(qset)
        total_feature_num = self.total_number_EmpiricalCompounds
        
        print("Query number of significant compounds = %d compounds" %query_set_size)
        
        for P in self.pathways:
            # use the measured pathway size
            P.overlap_EmpiricalCompounds = P.overlap_features = qset.intersection(P.EmpiricalCompounds)

            P.overlap_size = overlap_size = len(P.overlap_EmpiricalCompounds)
            P.EmpSize = ecpd_num = len(P.EmpiricalCompounds)
            if overlap_size > 0:
                negneg = total_feature_num + overlap_size - ecpd_num - query_set_size
                # Fisher's exact test
                P.p_FET = stats.fisher_exact([[overlap_size, query_set_size - overlap_size],
                                   [ecpd_num - overlap_size, negneg]], 'greater')[1]
                # EASE score as in Hosack et al 2003
                # taking out EASE, as the new approach of EmpiricalCompound is quite stringent already
                P.p_EASE = P.p_FET
                

            else:
                P.p_FET = P.p_EASE = 1
                
            FET_tested_pathways.append(P)
            #  (enrich_pvalue, overlap_size, overlap_features, P) 
            
        result = self.get_adjust_p_by_permutations(FET_tested_pathways)
        result.sort(key=lambda x: x.adjusted_p, reverse=False)
        self.resultListOfPathways = result

    
    def collect_hit_Trios(self):
        '''
        get [(mzFeature, EmpiricalCompound, cpd),...] for sig pathways.
        Nominate top cpd for EmpCpd here, i.e.
        in an EmpCpd promoted by a significant massFeature, the cpd candidate is chosen from a significant pathway.
        If more than one cpds are chosen, keep multiple.


        is this used now?
        main.AN = ActivityNetwork( mixedNetwork, set(PA.collect_hit_Trios() + MA.collect_hit_Trios()) )



        
        '''
        overlap_EmpiricalCompounds = set([])
        for P in self.resultListOfPathways:
            if P.adjusted_p < SIGNIFICANCE_CUTOFF:
                # print(P.adjusted_p, P.name)
                overlap_EmpiricalCompounds = overlap_EmpiricalCompounds.union(P.overlap_EmpiricalCompounds)
        
        new = []
        for T in self.TrioList:
            # [(mzFeature, EmpiricalCompound, cpd),...]
            if T[1] in overlap_EmpiricalCompounds and T[0] in self.mixedNetwork.significant_features:
                # this does not apply to all sig EmpCpd
                T[1].update_chosen_cpds(T[2])
                T[1].designate_face_cpd()
                new.append(T)
        
        return new
                    
    def to_json(self):
        '''
        Convert result to dataframes, easy JSON export to be consumed by downstream functions

        dict from cpd to empCpd:
        self.Compounds_to_EmpiricalCompounds - needs to clean up for pathway specific
        # 'dicts_cpd2empCpd': self.Compounds_to_EmpiricalCompounds,

        dicts_cpd2empCpd = []
        for P in self.resultListOfPathways:
            for E in P.overlap_EmpiricalCompounds:
                dicts_cpd2empCpd.append ({ E.chosen_compounds: E.EID })

        '''
        self.collect_hit_Trios()    # force update of overlap_EmpiricalCompounds

        L = []
        for P in self.resultListOfPathways:
            L.append({
                'pathway_id': P.id,
                'name': P.name, #.encode("utf-8", "ignore"),   # Pathway name not ascii compliant
                'overlap_size': P.overlap_size,
                'pathway_size': P.EmpSize,
                'p-value': P.adjusted_p ,
                'significant_empCpds': [ E.EID for E in P.overlap_EmpiricalCompounds],
                'significant_compounds': [";".join(E.chosen_compounds) for E in P.overlap_EmpiricalCompounds],
            })
        return L
