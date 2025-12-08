'''

    
'''

import json
import os
import logging
MASS_RANGE = (50, 2000)
RETENTION_TIME_TOLERANCE_FRAC = 0.02    


class InputUserData:
    '''
    
    Need to field 
    JSON list of features and list of epds
    
    User annotation if provided

    '''
    
    def __init__(self, paradict, web=False):
        '''
        
        '''
        self.web = web
        self.paradict = paradict
        self.header_fields = []
        self.ListOfUserFeatures = []
        self.input_featurelist = []

        # entry point of data input
        self.read()
        self.update()
        
    def update(self):
        '''
        Update retention_time_rank and is_significant to all MassFeatures
        
        Get user supplied JSON annotation if provided.
        
        '''
        retention_times = [M['rtime'] for M in self.ListOfUserFeatures]
        self.max_retention_time = max(retention_times)

        self.max_mz = max([M['mz'] for M in self.ListOfUserFeatures])
        self.determine_significant_list(self.ListOfUserFeatures)
        
        if 'annotation' in self.paradict and self.paradict['annotation']:
            # load empirical compound annotation from JSON file
            empCpd_json_file = os.path.join(self.paradict['workdir'], self.paradict['annotation'])
            with open(empCpd_json_file) as f:
                empCpd_json_text = f.read()
            self.EmpiricalCompounds = json.loads(empCpd_json_text)
            
            print("Loaded %d empirical compounds from annotation file." %len(self.EmpiricalCompounds))
            
        else:
            self.EmpiricalCompounds = {}
        
        
    def text_to_ListOfUserFeatures(self, textValue, delimiter='\t'):
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
                self.ListOfUserFeatures.append( 
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
        Read input feature lists to ListOfUserFeatures. 
        Row_numbers (rowii+1) are used as primary ID.
        # not using readlines() to avoid problem in processing some Mac files
        '''
        if self.web:
            self.text_to_ListOfUserFeatures(self.paradict['datatext'])
        else:
            self.text_to_ListOfUserFeatures( 
                open(os.path.join(self.paradict['workdir'], self.paradict['infile'])).read() )

        print("Read %d features as reference list." %len(self.ListOfUserFeatures))
    
    
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
        for f in self.ListOfUserFeatures:
            if f['pval'] < self.paradict['cutoff']:
                f['is_significant'] = True
            else:
                f['is_significant'] = False
        
        self.input_featurelist = [f['fid_from_user'] for f in self.ListOfUserFeatures if f['is_significant']]
        print("Using %d features (p < %f) as significant list." 
                              %(len(self.input_featurelist), self.paradict['cutoff']))  

