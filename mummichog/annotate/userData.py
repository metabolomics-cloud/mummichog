'''

0) JMS converts a JSON metabolic model to metabolicModels

1) JMS treats a metabolic model as a database, via KCD
2) JMS converts a list of features to empCpds in EED
3) return EED.dict_empCpds


from jms.modelConvert import convert_json_model
from jms.empiricalCpds import load_epds_from_json


# jms-metabolite-services
from jms.dbStructures import knownCompoundDatabase, ExperimentalEcpdDatabase

from .parameters import adduct_search_patterns, \
                                adduct_search_patterns_neg, \
                                    isotope_search_patterns, \
                                        extended_adducts


        EED = ExperimentalEcpdDatabase(mode=self.mode, mz_tolerance_ppm=self.mz_tolerance_ppm, rt_tolerance=2)
        # passing patterns from .defaul_parameters
        if self.mode == 'pos':
            EED.adduct_patterns = adduct_search_patterns
        else:
            EED.adduct_patterns = adduct_search_patterns_neg
        EED.isotope_search_patterns = isotope_search_patterns
        EED.extended_adducts = extended_adducts

        EED.build_from_list_peaks(self.CMAP.FeatureList)
        EED.extend_empCpd_annotation(self.KCD)
        EED.annotate_singletons(self.KCD)       
        # EED.dict_empCpds misses some features 
        EED.dict_empCpds = self.append_orphans_to_epmCpds(EED.dict_empCpds)


 [{'feature': {'rtime': 25.07, 'mz': 85.0647, 'id': 'F45'},
   'lib': [{'id': 'row199',
     'name': 'Cyclopentanone',
     'mz': 85.064791478,
     'rtime': 26.4}],
   'ms2': {},
   'csm': {}},
  {'feature': {'rtime': 191.32, 'mz': 85.0884, 'id': 'F65'},
   'lib': [],
   'ms2': {},
   'csm': {'name': 'piperidino',
    'mz': 85.0884,
    'rtime': 191.32,
    'ion_relation': nan,
    'empCpd_id': nan}},
  {'feature': {'rtime': 24.02, 'mz': 85.101, 'id': 'F66'},
   'lib': [],
   'ms2': {},
   'csm': {'name': 'Cyclohexane',
    'mz': 85.101,
    'rtime': 24.02,
    'ion_relation': nan,
    'empCpd_id': nan}}]
    
    
'''

import json
import os
import logging
MASS_RANGE = (50, 2000)
RETENTION_TIME_TOLERANCE_FRAC = 0.02    
# from mass2chem.adducts import *


#
# this should be based on Experiment class, and organize by peaks - features - empCpds
#

class InputUserData:
    '''
    
    Changing in v3
    need to field 
    JSON list of features and list of epds


    '''
    
    def __init__(self, paradict, web=False):
        '''
        
        '''
        self.web = web
        self.paradict = paradict
        self.header_fields = []
        self.ListOfMassFeatures = [] # MassFeatures are legacy name; they are just features
        self.input_featurelist = []

        # entry point of data input
        self.read()
        self.update()
        
    def update(self):
        '''
        Update retention_time_rank and is_significant to all MassFeatures
        '''
        retention_times = [M['rtime'] for M in self.ListOfMassFeatures]
        self.max_retention_time = max(retention_times)

        self.max_mz = max([M['mz'] for M in self.ListOfMassFeatures])
        self.determine_significant_list(self.ListOfMassFeatures)
        
        if 'annotation' in self.paradict and self.paradict['annotation']:
            # load empirical compound annotation from JSON file
            empCpd_json_file = os.path.join(self.paradict['workdir'], self.paradict['annotation'])
            with open(empCpd_json_file) as f:
                empCpd_json_text = f.read()
            self.EmpiricalCompounds = json.loads(empCpd_json_text)
            
            print("Loaded %d empirical compounds from annotation file." %len(self.EmpiricalCompounds))
            
        else:
            self.EmpiricalCompounds = {}
        
        
    def text_to_ListOfMassFeatures(self, textValue, delimiter='\t'):
        '''
        Column order is hard coded for now, as mz, retention_time, p_value, statistic, CompoundID_from_user

        use asari style JSON features

        '''
        def _make_id(ii, mz, rt):
            return 'F' + str(ii) + '_' + str(round(mz, 6)) + '@' + str(round(rt, 2))
        #
        lines = self.__check_redundant__( textValue.splitlines() )
        self.header_fields = lines[0].rstrip().split(delimiter)

        excluded_list = []
        for ii in range(len( lines )-1):
            y = lines[ii+1].split('\t')
            
            fid_from_user = ''
            if len(y) > 4: fid_from_user = y[4].strip()

            [mz, rtime, p_value, statistic] = [float(x) for x in y[:4]]
            
            # row_number, mz, retention_time, p_value, statistic, CompoundID_from_user
            if MASS_RANGE[0] < mz < MASS_RANGE[1]:
                # row # human-friendly, numbering from 1
                fid = _make_id(ii+1, mz, rtime)
                peak = {'id_number': fid, 
                        'id': fid,
                        'fid_from_user': fid_from_user or fid,
                        'mz': mz, 
                        'rtime': rtime,
                        'pval': p_value,
                        'statistic': statistic,
                        }
                self.ListOfMassFeatures.append( 
                    peak
                    )
            else:
                excluded_list.append( (ii, mz, rtime) )
        
        if excluded_list:
            print( "Excluding %d features out of m/z range %s." %(len(excluded_list), str(MASS_RANGE)) )

        
    def read_from_file(self, inputFile):
        return open(inputFile).read()
    
    def read_from_webform(self, t):
        return t

    def __check_redundant__(self, L):
        redundant = len(L) - len(set(L))
        if redundant > 0:
            print( "Your input file contains %d redundant features." %(redundant) )
        return L

    def read(self):
        '''
        Read input feature lists to ListOfMassFeatures. 
        Row_numbers (rowii+1) are used as primary ID.
        # not using readlines() to avoid problem in processing some Mac files
        '''
        if self.web:
            self.text_to_ListOfMassFeatures(self.paradict['datatext'])
        else:
            self.text_to_ListOfMassFeatures( 
                open(os.path.join(self.paradict['workdir'], self.paradict['infile'])).read() )

        print("Read %d features as reference list." %len(self.ListOfMassFeatures))
    
    
    # more work?
    def determine_significant_list(self, all_feature_list):
        '''
        For single input file format in ver 2. 
        The significant list, input_mzlist, should be a subset of ref_mzlist,
        determined either by user specificed --cutoff,
        or by automated cutoff close to a p-value hotspot,
        in which case, paradict['cutoff'] is updated accordingly.

        '''
        if not self.paradict['cutoff']:
            # automated cutoff
            new = sorted(all_feature_list, key=lambda x: x['pval'])
            
            p_hotspots = [ 0.2, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0001 ]
            N_hotspots = [ len([x for x in all_feature_list if x['pval'] < pp]) for pp in p_hotspots ]
            
            N_quantile = len(new) / 4
            N_optimum, N_minimum = 300, 30
            chosen = 9999
            for ii in range( len(N_hotspots) ):
                # will get the smallest p as ii increases
                if N_optimum < N_hotspots[ii] < N_quantile:
                    chosen = ii
            
            # if nothing was chosen
            if chosen > 100:
                for ii in range( len(N_hotspots) ):
                    if N_minimum < N_hotspots[ii] < N_quantile:
                        chosen = ii
            
            if chosen > 100:
                N_chosen = int(N_quantile)
                self.paradict['cutoff'] = new[N_chosen+1]['pval']
            else:
                #N_chosen = N_hotspots[chosen]
                
                self.paradict['cutoff'] = p_hotspots[chosen]
        
            print("Automatically choosing (p < %f) as significant cutoff."  %self.paradict['cutoff'])  
        
        # mark MassFeature significant
        for f in self.ListOfMassFeatures:
            if f['pval'] < self.paradict['cutoff']:
                f['is_significant'] = True
            else:
                f['is_significant'] = False
        
        self.input_featurelist = [f['fid_from_user'] for f in self.ListOfMassFeatures if f['is_significant']]
        print("Using %d features (p < %f) as significant list." 
                              %(len(self.input_featurelist), self.paradict['cutoff']))  

