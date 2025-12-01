'''
Example of empirical compound JSON format:

{
"interim_id": "kp1832_268.0808",
"neutral_formula_mass": 268.08077, 
"neutral_formula": "C10H12N4O5",
"identity": [
        {"compounds": ["HMDB0000195"], "names": ["Inosine"], 
                "score": 0.6, "probability": null},
        {"compounds": ["HMDB0000195", "HMDB0000481"], "names": ["Inosine", "Allopurinol riboside"], 
                "score": 0.1, "probability": null},
        {"compounds": ["HMDB0000481"], "names": ["Allopurinol riboside"], 
                "score": 0.1, "probability": null},
        {"compounds": ["HMDB0003040"], "names": ["Arabinosylhypoxanthine"], 
                "score": 0.05, "probability": null}
        ],
"MS1_pseudo_Spectra": [
        {"feature_id": "FT1705", "mz": 269.0878, "rtime": 99.90, 
                "isotope": "M0", "modification": "+H", "ion_relation": "M+H[1+]"},
        {"feature_id": "FT1876", "mz": 291.0697, "rtime": 99.53, 
                "isotope": "M0", "modification": "+Na", "ion_relation": "M+Na[1+]"},
        {"feature_id": "FT1721", "mz": 270.0912, "rtime": 99.91, 
                "isotope": "13C", "modification": "+H", "ion_relation": "M(C13)+H[1+]"},
        {"feature_id": "FT1993", "mz": 307.0436, "rtime": 99.79, 
                "isotope": "M0", "modification": "+K", "ion_relation": "M+K[1+]"},
        ],
"MS2_Spectra": [
        {"precursor_ion_mz": 269.0879, "matched_MS1_feature": "FT1705", 
        "precursor_ion_id": "269.087921142578_100.9166816170002_plasma_ID_01.mzML",
        "source": "plasma_ID_01.mzML"},
        {"precursor_ion_mz": 307.0438, "matched_MS1_feature": "FT1993", 
        "precursor_ion_id": "307.043762207031_99.3396262240002_plasma_ID_05.mzML",
        "source": "plasma_ID_05.mzML"},
        ],
"annotation": {
        "authLib_Li_Lab_202410": [], 
        "MoNA-202402": [{
            "precursor_ion_id": "269.087921142578_100.9166816170002_plasma_ID_01.mzML",
            "reference_id": "Inosine",  "db_precursor_mz": 269.088, "db_precursor_rt": 43.44,
            "msms_score": 0.985}
        ],
        "HMDBv5": [{"accession": "HMDB0000195", "name": "Inosine"},
            {"accession": "HMDB0000481", "name": "Allopurinol riboside"},
            {"accession": "HMDB0003040", "name": "Arabinosylhypoxanthine"},
            {"accession": "HMDB0252449", "name": "Formycin b"}], 
        "CSM_r1": []   
        },
"Databases_referred": ["CSM_r1", "HMDBv5", "MoNA-202402", "authLib_Li_Lab_202410"]
}


'''
import json

MASS_RANGE = (50, 2000)
RETENTION_TIME_TOLERANCE_FRAC = 0.02    


def score_cpd_identity(empCpd):
    '''
    Given identity list from empirical compound, score each compound based on 
    its identity entries.
    
    Equal weight for now. 
    CSM has a function (meant for factor graph) that can be used here.
    
    We need cpd ID standardization here;
    then make sure they are consistent with metabolic model.
        
    Return:
        {cpd_id: score, ...}
    
    '''
    cpd_scores = {}
    compounds = []
    # Below is a hack; need proper cpd count & scoring function later
    if empCpd['identity']:
        for x in empCpd['identity']:
            if len(x['compounds']) == 1:    # ignore mixtures for now
                cpd_scores[x['compounds'][0]] = x.get('score', 0.1)
    # need ID standardization here
    else:
        for k,v in empCpd.get('annotation', {}).items():
            if k.startswith('HMDB'):
                compounds += [entry.get('accession') for entry in v]
            elif k.startswith('KEGG'):
                compounds += [entry.get('accession') for entry in v]
            elif k.startswith('authLib_Li_Lab'):
                compounds += [entry.get('cpd') for entry in v]
            elif k.startswith('MoNA'):
                compounds += [entry.get('reference_id') for entry in v]
            # add more databases here as needed
        for cpd in set(compounds):
            if cpd not in cpd_scores:
                cpd_scores[cpd] = 0.1    # default score for unscored IDs

    return cpd_scores


class DataMeetModel:
    '''
    TrioList format can be still used in v3

    '''
        
    def __init__(self, metabolicModel, userData):
        '''
        # from ver 1 to ver 2, major change in .match()
        Trio structure of mapping
        (FeatureID, EmpiricalCompounds, Cpd)
        
        '''
        self.model = metabolicModel
        self.data = userData
        self.rowDict = {f['id']: f for f in self.data.ListOfMassFeatures}
        # this is the sig list
        self.significant_features = self.data.input_featurelist
        # this is the reference list
        self.features = [f['id'] for f in self.data.ListOfMassFeatures]
        self.ListOfEmpiricalCompounds = self.List_score_EmpiricalCompounds()
        
        self.feature_to_EmpiricalCompound, self.Compound_to_EmpiricalCompounds, \
            self.TrioList = self.index_EmpCpd_Cpd()

            
    def List_score_EmpiricalCompounds(self):
        '''
        EmpiricalCompounds should be already constructed in userData. 
        If not, do it here.
        
        Singletons may have matched neutral_formula using primary ions. 
        If not, deal wtih singletons, as M+H+ or M-H- forms.
        
        Need to make sure cpd IDs are consistent with those in metabolic model.
        '''
        
        
        if not self.data.EmpiricalCompounds:
            # construct from scratch using khipu
            
            pass
        
        for empCpd in self.data.EmpiricalCompounds.values():
            # first deal with singletons
            if not empCpd['neutral_formula_mass']:
                empCpd['cpd_scores'] = {}
                mz = empCpd['MS1_pseudo_Spectra'][0]['mz']
                if self.data.paradict.get('ionization', 'pos') == 'pos':
                    empCpd['neutral_formula_mass'] = mz - 1.007276
                else:
                    empCpd['neutral_formula_mass'] = mz + 1.007276
                
            else: # khipu should have some annotation already
                # {cpd_id: score, ...}
                empCpd['cpd_scores'] = score_cpd_identity(empCpd)
                
            # add new round of matching to metabolic model here
            self.augment_empCpd_with_model_cpds(empCpd)
        
        return self.data.EmpiricalCompounds


    def augment_empCpd_with_model_cpds(self, empCpd):
        '''
        Given an empirical compound, augment its identity list with compounds
        from metabolic model, based on neutral_formula_mass matching.
        
        Update empCpd in place.
        
        '''
        ppm = self.data.paradict['ppm'] or 10
        mass_tol = ppm * empCpd['neutral_formula_mass'] / 1e6
        matched_cpds = []
        #
        # yet to sort this out
        #
        for cpd_id, cpd in self.model.Compounds.items():
            if abs(cpd.neutral_formula_mass - empCpd['neutral_formula_mass']) <= mass_tol:
                matched_cpds.append(cpd_id)
        
        # add matched_cpds to identity list with default score
        for cpd_id in matched_cpds:
            if cpd_id not in empCpd['cpd_scores']:
                empCpd['cpd_scores'][cpd_id] = 0.05    # default score for model-matched IDs


    def index_EmpCpd_Cpd(self):
        '''
        mapping from feature row index to EmpiricalCompounds and Compounds.
        
        self.ListOfEmpiricalCompounds may contain more features than self.features.
        Only those in self.features are used to build TrioList.
        
        Cpd IDs have to be consistent with those in metabolic model (self.get_ListOfEmpiricalCompounds).
        Singletons are matched to model as M+H+ or M-H-.
        
        Return:
            feature_to_EmpiricalCompounds: {feature: [EmpiricalCompounds, ...], ...}
            Compounds_to_EmpiricalCompounds: {cpd: [EmpiricalCompounds, ...], ...}
            TrioList: [(feature, EmpiricalCompounds, cpd), ...]
        
        '''
        feature_to_EmpiricalCompounds, cpd2EmpiricalCompounds = {}, {}
        TrioList = []
        for empCpd in self.ListOfEmpiricalCompounds:
            for feat in empCpd['MS1_pseudo_Spectra']:
                feature_to_EmpiricalCompounds[feat['id']] = empCpd['interim_id']
            for cpd, score in empCpd['cpd_scores'].items():
                if cpd not in cpd2EmpiricalCompounds:
                    cpd2EmpiricalCompounds[cpd] = []
                cpd2EmpiricalCompounds[cpd] += [empCpd['interim_id']]
                for feat in empCpd['MS1_pseudo_Spectra']:
                    TrioList.append( (feat['id'], empCpd['interim_id'], cpd) )
                    
        return feature_to_EmpiricalCompounds, cpd2EmpiricalCompounds, TrioList

    
    def to_json(self):
        '''
        JSON export to be consumed by downstream functions

        empCpd2Cpds = {empCpd: (), ...,}

        Will update later in accordance to 
        https://github.com/shuzhao-li/metDataModel

        '''

        empCpd2Features, empCpd2Cpds = {}, {}
        for E in self.ListOfEmpiricalCompounds:
            empCpd2Features[E.EID] = E.massfeature_rows
            empCpd2Cpds[E.EID] = E.compounds

        return {
            'metabolic_model': self.model.version,
            'empCpd2Features': empCpd2Features,
            'empCpd2Cpds': empCpd2Cpds,
        }
