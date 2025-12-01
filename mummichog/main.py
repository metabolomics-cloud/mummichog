# Licensed under the BSD 3-Clause License.
#
# mummichog - pathway and network analysis for metabolomics
#

import time
import argparse
import sys
import json
from mummichog import __version__
from mummichog.models.get_models import get_metabolic_model

from .api import *
from .parameters import PARAMETERS

fishlogo = '''     
    --------------------------------------------
    
             oO                      ooooooooo
           oOO   OOOOO  ooooo       ooo oooo
     oOO   O       ooooo  oooooo ooooo
    oooO           oooooo         oooo ooooo
        Oooo   o      OOOOOO   oooo   oooooooo
            ooooo  oooo      
                 o
    
    --------------------------------------------
    '''

def build_parser():
    parser = argparse.ArgumentParser(
        description='mummichog v%s: pathway and network analysis for metabolomics' %__version__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # add arguments
    parser.add_argument('-v', '--version', action='version', version=__version__, 
            help='print version and exit')
    parser.add_argument('-j', '--project', type=str,
            help='project name')
    parser.add_argument('-m', '--mode', type=str,
            help='mode of ionization, pos or neg')
    parser.add_argument('--ppm', type=int, 
            help='mass precision in ppm (part per million), same as mz_tolerance_ppm')
    parser.add_argument('-d', '--workdir', type=str,
            help='working directory')
    parser.add_argument('-i', '--infile', type=str,
            help='input file with statistical results')
    parser.add_argument('-a', '--annotation', type=str,
            help='annotation file in empirical compound format (json)')
    parser.add_argument('-o', '--output', type=str,
            help='output directory')
    parser.add_argument('-c', '--cutoff', type=float,
            help='significance cutoff for features, p-value or similar metric')

    parser.add_argument('-p', '--permutation', type=int,
            help='number of permutations to estimate null distributions')
    
    args = parser.parse_args()
    return args


def main():
    
    print (fishlogo)
    print ( "mummichog version %s \n" %__version__ )
    
    # make a copy of the default parameters; user options will override
    parameters = PARAMETERS.copy()

    # build CLI parser
    args = build_parser()
    parameters.update(vars(args)) 

    print("Started @ %s\n" %time.asctime())
    userData = InputUserData(parameters)
    theoreticalModel = get_metabolic_model( parameters['network'] )
    
    # for developer testing
    print(
            list(theoreticalModel.Compounds.items())[92], "...\n"
    )
    print(parameters)
    
    mixedNetwork = DataMeetModel(theoreticalModel, userData)
    
    
    # getting a list of Pathway instances, with p-values, in PA.resultListOfPathways
    PA = PathwayAnalysis(mixedNetwork.model.metabolic_pathways, mixedNetwork)
    PA.cpd_enrich_test()
    
    # Module analysis, getting a list of Mmodule instances
    MA = ModularAnalysis(mixedNetwork)
    MA.dispatch()
    
    # do activity network
    AN = ActivityNetwork( mixedNetwork, set(PA.collect_hit_Trios() + MA.collect_hit_Trios()) )


    print("\nFinished @ %s\n" %time.asctime())

    #
    #  This is to export data as Python objects
    #
    MCG_JSON = json_export_all(mixedNetwork, PA, MA, AN)

    #print(MCG_JSON)
    print("\n\n~~~~~~~~~~~~~~~~~~~~\n\n")

    # Use Encoder to convert Python to JSON format 
    s = json.JSONEncoder().encode(MCG_JSON )
    with open("mcg_output.json", "w") as O:
        O.write(s)

    #print(json.dumps(s,  indent=4) [:500])
    print("JSON output was written in mcg_output.json.")

  

#
# -----------------------------------------------------------------------------
#

if __name__ == '__main__':

    main()


