'''Mummichog parameter templates
in progress
Pre-annotation will be moved out of mummichog, but based on khipu.

'''

PROTON = 1.00727646677
electron = 0.000549

SIGNIFICANCE_CUTOFF = 0.05
MASS_RANGE = (50, 2000)

# fraction of total retention time, or of ranks of retention time
# used to determine coelution of ions ad hoc
RETENTION_TIME_TOLERANCE_FRAC = 0.02    

PARAMETERS = {
    'network': 'human_model_mfn',   # metabolic model to use
    'mode': 'pos_default',    # analytical mode of mass spec
    'instrument': 'unspecified',  # instrument type, for future use
    'cutoff': 0.05,              # p-value cutoff to select significant features
    'ppm': 10,                  # mass precision in ppm (part per million)
    'modeling': None,         # modeling permutation data, None or 'gamma'
    'input': '',              # input data file
    'output': '',             # output file prefix
    'permutation': 100,       # number of permutations to estimate null distributions
    'outdir': 'mcgresult',    # output directory name
}
