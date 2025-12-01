'''
import csv
import xlsxwriter
import logging

'''

def json_export_all(mixedNetwork, PA, MA, AN):
    '''
    metabolic model is already in JSON, but need clean up.
    '''
    return {
        'EmpiricalCompounds': mixedNetwork.to_json(),
        'pathway_analysis': PA.to_json(),     #force_ascii=True),
        'module_analysis': MA.to_json(), 
        'activity_network': AN.to_json(),
    }
    


# Not used below for now

# --------------------------------------------------------
#
# plot in pathway analysis, self=PA
#


    def plot_model_pvalues(self, outfile='mcg_pathway_modeling'):
        '''
        Plot self.permutation_record
        P.p_EASE for P in self.resultListOfPathways
        Use -log10 scale, to show upward trend, consistent with other plots
        '''
        self.permutation_record.sort()
        Y_data = [-np.log10(x) for x in self.permutation_record]
        fig = plt.figure(figsize=(5,4))
        plt.plot(range(len(Y_data)), Y_data, 'b.')
        for P in self.resultListOfPathways[:10]:
            YY = -np.log10(P.p_EASE)
            plt.plot([0, 0.1*len(Y_data)], [YY, YY], 'r-')
        
        plt.ylabel("-log10 (FET p-value)")
        plt.xlabel("Number of permutation")
        plt.title("Modeling pathway significance")
        plt.tight_layout()
        plt.savefig(outfile+'.pdf')

    
    def plot_bars_top_pathways(self, outfile='mcg_pathway_barplot'):
        '''
        Horizontal barplot of pathways.
        Also returnin-memory string for web use
        '''
        use_pathways = [P for P in self.resultListOfPathways if P.adjusted_p < SIGNIFICANCE_CUTOFF]
        if len(use_pathways) < 6:
            use_pathways = self.resultListOfPathways[:6]
        #plot use_pathways
        fig, ax = plt.subplots()
        ylabels = [P.name for P in use_pathways]
        data = [-np.log10(P.adjusted_p) for P in use_pathways]
        NN = len(data)
        ax.barh( range(NN), data, height=0.5, align='center', color="purple", alpha=0.4 )
        ax.set_yticks(range(NN))
        ax.set_yticklabels(ylabels)
        ax.set_xlabel('-log10 p-value')
        
        ax.plot([1.301, 1.301], [-0.5, NN], 'g--')    # NN is inverted too
        #ax.set_ylim(-0.5, NN)
        
        ax.invert_yaxis()
        plt.tight_layout()
        plt.savefig(outfile+'.pdf')
        
        # get in-memory string for web use
        figdata = BytesIO()
        plt.savefig(figdata, format='png')
        return """<img src="data:image/png;base64,{}"/>""".format(base64.encodebytes(figdata.getvalue()).decode()) 


# --------------------------------------------------------
#
# plot in module analysis, self=MA
#

    def plot_model_pvalues(self, outfile='mcg_module_modeling.pdf'):
        '''
        Plot module activity against self.permuation_mscores
        
        '''
        self.permuation_mscores.sort(reverse=True)
        NN = len(self.permuation_mscores)
        fig = plt.figure(figsize=(5,4))
        plt.plot(range(NN), self.permuation_mscores, 'bo')
        for M in self.modules_from_significant_features:
            plt.plot([0, 0.1*NN], [M.A, M.A], 'r-')
        
        plt.ylabel("Activity score")
        plt.xlabel("Number of permutation")
        plt.title("Modeling module significance")
        plt.tight_layout()
        plt.savefig(outfile+'.pdf')
        
    
    def draw_top_modules(self, outfile_prefix='mcg_module_.pdf'):
        for M in self.top_modules:
            #draw it
            #
            pass
    