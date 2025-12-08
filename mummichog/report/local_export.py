'''
Old export functions

Not needed in v3 web app,
but can be optionally used in local and in Jupyter Notebook applications


'''

import os
import csv

from .websnippets import *

#
# to remove dependency from matplotlib in ver 3
#
import matplotlib.pyplot as plt

class WebReporting:
    '''
    To generate HTML summary report.
    
    This will change to template (Jinja2) based methods.
    '''
    def __init__(self, Local, PA, MA, AN):

        self.mixedNetwork = Local.mixedNetwork
        self.model = Local.mixedNetwork.model
        self.data = Local.mixedNetwork.data
        self.Local = Local
        
        self.PA = PA
        self.MA = MA
        self.AN = AN
        
    def run(self):
        self.get_dict_cpd_statistic()
        self.collect_web_export_graphs()
        html = self.web_export()
        return html
        
    def run_nodisk(self):
        '''
        For web use, no local disk write
        '''
        self.get_dict_cpd_statistic()
        self.collect_web_export_graphs()
        return self.web_export_str()
        
    def get_dict_cpd_statistic(self):
        self.dict_cpd_statistic = {}
        for T in self.mixedNetwork.hit_Trios:
            self.dict_cpd_statistic[T[2]] = self.mixedNetwork.rowDict[T[0]].statistic
        
    def collect_web_export_graphs(self):
        '''
        Return activity network and up to top 5 modules.
        Moving activity network after modules
        
        AN yet to be revised
        
        '''
        self.web_export_networks = [M.graph for M in self.MA.top_modules[:5]] + [self.AN.activity_network]
        
        
        vis_nodes = []
        for g in self.web_export_networks: vis_nodes += g.nodes()
        self.vis_nodes = set(vis_nodes)


    def web_export(self):
        '''
        Return HTML without head/foot.
        '''
        outfile = os.path.join(self.Local.rootdir, 'result.html')
        webtextList = self.web_export_str()
        with open(outfile, 'w') as O:
            O.write(''.join(webtextList))
        return webtextList[1]

    def web_export_str(self):
        '''
        Write HTML and javascript based report.
        The visualization section includes activity network and up to top 5 modules.
        return webtextList as [head, middle HTML, end]
        '''

        HTML = HtmlExport()
        title = 'Mummichog Report: ' + self.data.paradict['output']
        HTML.add_element(title, 'h1', '')
        
        HTML.add_element('', 'div', 'colorbar', 'colorbar')
        stats = "mummichog version: %s, using metabolic model %s.<br>" %(VERSION, self.model.version)
        HTML.add_element(stats, 'div', 'stats')
        
        # MWAS figure
        HTML.add_element('User input data', 'h2', '')
        details_userData = "User supplied %d features as reference list." %len(self.data.ListOfMassFeatures)
        details_userData += "We are using %d features (p < %f) as significant list. The feature level data are shown in the Manhattan plots below." %(
                                    len(self.data.input_featurelist), self.data.paradict['cutoff'] )
        HTML.add_element(details_userData, 'p', '')
        HTML.add_element(self.Local.inline_plot_userData_MWAS, 'div', 'inline_plot_userData_MWAS')

        # Pathway table and figure
        HTML.add_element('Top pathways', 'h2', '')
        details_pathwayData = "The reference feature list is mapped to %d EmpiricalCompounds, significant features to %d EmpiricalCompounds. " %(
                                    self.PA.total_number_EmpiricalCompounds, len(self.PA.significant_EmpiricalCompounds))
        details_pathwayData += "In the table below, pathway_size is number of detected EmpiricalCompounds for each pathway; overlap_size is number of significant EmpiricalCompounds.\
                                Empirical p-values are estimated by permutation test. Details on EmpiricalCompounds are in 'tables/ListOfEmpiricalCompounds.tsv'."
        HTML.add_element(details_pathwayData, 'p', '')
        pathwaystablly = self.write_pathway_table()
        HTML.add_element(pathwaystablly, 'div', 'pathwaystablly')
        HTML.add_element(self.Local.inline_plot_pathwayBars, 'div', 'inline_plot_pathwayBars')

        # place to insert network visusalization
        HTML.add_element('Top modules', 'h2', '')
        Module_pharase = "Network modules may capture concerted metabolite activities that are missed by predefined pathways. \
        Use pull-down menus to visualize modules and the activity network, which tries to combine the results of \
        pathway/module analyses. The visualization requires interenet connection."
        HTML.add_element(Module_pharase, 'p', '')
        select_menu = HTML.make_select_menu(len(self.web_export_networks))
        HTML.add_element(select_menu, 'div', 'network_selection', 'network_selection')
        # another color bar?
        #HTML.add_element('', 'div', 'colorbar', 'colorbar')
        HTML.add_element('', 'div', 'networkvisual', 'networkvisual')
        

        modulestablly = self.write_module_table()
        HTML.add_element(modulestablly, 'div', 'modulestablly')
        
        HTML.add_element('EmpiricalCompounds favored by above analyses', 'h2', '')
        EmpCpd_phrase = "An EmpiricalCompound is a computational unit for a tentative metabolite. It can group multiple ions, \
                    and be any of the isobaric/isomeric species. This table contains EmpiricalCompounds prioritized by pathway/module analysis."
        HTML.add_element(EmpCpd_phrase, 'p', '')
        metabolitestablly = self.write_metabolite_table()
        HTML.add_element(metabolitestablly, 'div', 'metabolitestablly')
        
        HTML.add_element('More data', 'h2', '')
        moreinformation = '''Full data tables from pathway analysis, module analysis and activity network are stored in the 'tables/' directory. 
            The tsv files can be imported into a spreadsheet program, e.g. MS Excel.
            Under the 'figures/network_modules/' directory are files intended for Cytoscape (cytoscape.org) visualization. Please refer to Cytoscape's guides for details.
            Details of this run are recorded in file mummichog.log.txt.
            '''
        HTML.add_element(moreinformation, 'p', 'moreinformation')
        
        footer = '''Mummichog algorithms are described in Li et al. Predicting Network Activity from High Throughput Metabolomics. PLoS Computational Biology (2013); doi:10.1371/journal.pcbi.1003123.\
        This software is provided as is, without warranty of any kind. Visit <a href="http://mummichog.org">http://mummichog.org</a> for updates.
        '''
        HTML.add_element(footer, 'footer', '')
        
        # moved to self.export_sif_related()
        # self.export_cpd_attributes()
        
        # push network data into javascript
        HTML.make_js_data(self.web_export_networks, 
                          self.model.dict_cpds_def, 
                          self.dict_cpd_statistic)
        
        return HTML.export_text()
        
    def filter_vis_nodesdict(self, d1):
        d2 = {}
        for k in d1.keys():
            if k in self.vis_nodes:
                d2[k] = d1[k]
                
        return d2

    def write_pathway_table(self):
        s = '<table><tr><th>Pathways</th>\
                    <th>overlap_size</th><th>pathway_size</th><th>p-value</th>\
                    <th>overlap_EmpiricalCompounds</th></tr>'
        ii = 0
        for P in self.PA.resultListOfPathways:
            ii += 1
            if ii < 6 or P.adjusted_p < SIGNIFICANCE_CUTOFF:
                empCpds = [E.EID for E in P.overlap_EmpiricalCompounds]
                s += '<tr> <td>' + P.name + '</td> <td>'\
                        + str(P.overlap_size) + '</td><td>' + str(P.EmpSize
                        ) + '</td><td>' + str(round(P.adjusted_p, 5)) + '</td><td>' \
                        + ','.join(empCpds) + '</td></tr>'        

        return s + '</table>'

    
    def write_module_table(self):
        '''
        for M in self.MA.top_modules:
            counter += 1
            M.export_network_txt(self.model, 
                                 os.path.join(self.moduledir, 'module_' + str(counter) + '.txt'))
        
        '''
        s, counter = '', 0
        for M in self.MA.top_modules:
            counter += 1
            nodes = M.graph.nodes()
            names = [self.model.dict_cpds_def.get(x, '') for x in nodes]
            s += '<div class="moduleline">' + 'module_' + str(counter) + ", p=" + str(round(M.p_value, 5)) + ", " + str(len(nodes)) + " metabolites" + '</div>'
            s += '<div class="metabolites">' + ', '.join(names) + '</div>'
            
        return s + '\n'



    def write_metabolite_table(self):
        '''
        To export significant EmpiricalCompounds that has a face_compound, which was added during pathway/module analysis
        
        
        '''
        
        s = '<table><tr><th>EmpiricalCompound</th><th>Input m/z</th>\
                    <th>Retention time</th><th>ion</th><th>mz_diff</th><th>Statistic</th>\
                    <th>Significant</th></tr>'
        # will add theoretical MW?
        for E in self.PA.significant_EmpiricalCompounds:
            if E.face_compound:
                col2 = 'Best guess: ' + E.face_compound + ', ' + self.model.dict_cpds_def.get(E.face_compound, '').split(";")[0]
                s += '<tr > <td>' + E.EID + '</td> <td colspan="2">' + col2 + '</td><td colspan="4">' + ', '.join(E.compounds) + '</td></tr>'
                for f in E.massfeature_rows:
                    F = self.mixedNetwork.rowDict[f]
                    ion = E.row_to_ion[f]
                    
                    # should store earlier?
                    mz_diff = F.mz - self.model.Compounds[ E.face_compound ]['adducts'][ ion ]
                    
                    s += '<tr> <td> </td> <td>' + str(F.mz) + '</td> <td>' + str(F.retention_time) + '</td><td>' + ion + '</td><td>' + str(
                                round(mz_diff,4)) + '</td><td>' + str(round(F.statistic,2)) + '</td><td>' + write_yes_no_MassFeature(F) + '</td></tr>'
           
        return s + '</table>'
                


class LocalExporting:
    '''
    To write out local data

    '''
    def __init__(self, mixedNetwork, PA, MA, AN):
        self.mixedNetwork = mixedNetwork
        self.model = mixedNetwork.model
        self.data = mixedNetwork.data
        
        self.PA = PA
        self.MA = MA
        self.AN = AN
        
        self.create_dirs()
        
    def create_dirs(self):
        '''
        Directory tree, renamed for version 2
        outdir/ result.html
            mummichog.log.txt
            tables/
            figures/
                network_modules/ 
        
        logging.basicConfig(filename=os.path.join(optdict['workdir'], optdict['outdir'], 'mummichog.log'), 
                            format='%(message)s', 
                            level=logging.INFO)
        logging.info('\n'.join(["mummichog version: %s" %VERSION,
                                 "pwd: %s" %os.getcwd(),
                                 "user command: %s" %' '.join(sys.argv),
                                 "\n",
                                 ]))
        
        
        #print("Pathway analysis report was written to \n%s (.tsv and .xlsx)" %outfile)
     


        Web export may avoid static file linking, as it's tricky to be correct on local computers
        '''
        self.rootdir = os.path.join(self.data.paradict['workdir'], self.data.paradict['outdir'])
        os.mkdir(self.rootdir)
        self.tabledir, self.figuredir, self.moduledir = os.path.join(
                                                self.rootdir, 'tables'), os.path.join(
                                                self.rootdir, 'figures'), os.path.join(
                                                self.rootdir, 'figures', 'network_modules')
        os.mkdir(self.tabledir)
        os.mkdir(self.figuredir)
        os.mkdir(self.moduledir)

    def run(self):
        '''
        '''
        # export tables
        self.export_userData()
        self.export_EmpiricalCompounds()
        self.export_pathway_enrichtest()
        self.writeTable_top_modules()
        
        # plot figures
        self.plot_userData_MWAS()
        self.plot_pathwayBars()
        self.plot_pathway_model()
        self.plot_module_model()
        
        # export mudules and activity network
        self.export_top_modules()
        self.export_activity_network()
        self.export_cpd_attributes()
        
        
    def export_userData(self):
        '''
        to do
        
        Should add this function to write out user data with row_numbers,
        to help users to track data.
        
        self.ListOfMassFeatures.append( 
                    MassFeature('row'+str(ii+1), mz, retention_time, p_value, statistic, CompoundID_from_user) )
        
        '''
        s = "massfeature_rows\tm/z\tretention_time\tp_value\tstatistic\tCompoundID_from_user\n"
        for F in self.data.ListOfMassFeatures:
            s += F.make_str_output() + '\n'
                
        with open(os.path.join(self.tabledir, "userInputData.txt"), 'w') as O:
            O.write(s)
            
            

    def export_EmpiricalCompounds(self):
        '''
        This exports all tri relationships.
        In ActivityNetwork, the top predicted metaolite is determined.

        To include the annotation model here, i.e. EmpiricalCompound.identity


        '''
        s = "EID\tmassfeature_rows\tstr_row_ion\tcompounds\tcompound_names\n"
        for E in self.mixedNetwork.DictOfEmpiricalCompounds.values():
            names = [self.model.dict_cpds_def.get(x, '') for x in E.compounds]
            s += '\t'.join([E.EID, ';'.join(E.massfeature_rows), E.str_row_ion, ';'.join(E.compounds), '$'.join(names)]
                ) + '\n'
        with open(os.path.join(self.tabledir, "ListOfEmpiricalCompounds.tsv"), 'w') as O:
            O.write(s)

    
    # export functions for pathway analysis
    def export_pathway_enrichtest(self):
        '''
        result as sorted list of pathway instances
        Cutting out 'p-value (raw)' from version 1.
        '''
        resultstr = [['pathway', 'overlap_size', 'pathway_size', 'p-value', 
                      'overlap_EmpiricalCompounds (id)', 'overlap_features (id)', 'overlap_features (name)',] ]
        
        for P in self.PA.resultListOfPathways:
            empCpds = [E.EID for E in P.overlap_EmpiricalCompounds]
            cpds = [E.chosen_compounds for E in P.overlap_EmpiricalCompounds]
            names = [ [self.model.dict_cpds_def.get(x, '') for x in y] for y in cpds ]
            resultstr.append([str(x) for x in [P.name, P.overlap_size, P.EmpSize, P.adjusted_p]]
                             + [','.join(empCpds), ','.join(['/'.join(x) for x in cpds]), '$'.join(['/'.join(x) for x in names]) ])

        outfile = os.path.join(self.tabledir, "mcg_pathwayanalysis_") + self.data.paradict['output']
        # write .tsv
        f =  open(outfile + '.tsv', "w", newline="")
        writer = csv.writer( f, delimiter='\t' )
        writer.writerows( resultstr )
        f.close()
        
        # export .xlsx too
        write_xlsx(resultstr, outfile + '.xlsx', "mcg_pathwayanalysis_")
        

    # export functions for module analysis
    def writeTable_top_modules(self):
        '''
        Write .tsv report for top modules.
        '''
        resultstr = [
                     ['MODULE', 'p-value', 'size', 'members (name)', 'members (id)', 
                      'This module overlaps with'], ]
        counter = 0
        for M in self.MA.top_modules:
            counter += 1
            nodes = M.graph.nodes()
            names = [self.model.dict_cpds_def.get(x, '') for x in nodes]
            resultstr.append(['module_' + str(counter), str(M.p_value),
                              str(len(nodes)), '$'.join(names),
                              ','.join(nodes), self.find_top_pathways(nodes),
                               ])

        outfile = os.path.join(self.tabledir, "mcg_modularanalysis_") + self.data.paradict['output']
        # write .tsv
        f =  open(outfile + '.tsv', "w", newline="")
        writer = csv.writer( f, delimiter='\t' )
        writer.writerows( resultstr )
        f.close()
        # export .xlsx too
        write_xlsx(resultstr, outfile + '.xlsx', "mcg_modularanalysis_")
        
    def find_top_pathways(self, nodes):
        '''
        return top 3 pathways for each feature, used for annotation
        '''
        pathways = []
        for n in nodes: pathways += self.model.cpd2pathways.get(n, [])
        counts = [(pathways.count(p), p) for p in set(pathways)]
        counts.sort(reverse=True)
        return ';'.join([str(x) for x in counts[:3]])


    def plot_userData_MWAS(self):
        self.inline_plot_userData_MWAS = self.data.make_manhattan_plots( 
            os.path.join(self.figuredir, "mcg_MWAS_") + self.data.paradict['output']
            )

    def plot_pathwayBars(self):
        self.inline_plot_pathwayBars = self.PA.plot_bars_top_pathways(
            os.path.join(self.figuredir, "mcg_pathwayBars_") + self.data.paradict['output']
            )

    def plot_pathway_model(self):
        self.PA.plot_model_pvalues(
            os.path.join(self.figuredir, "plot_pathwayModel_") + self.data.paradict['output']
            )

    def plot_module_model(self):
        self.MA.plot_model_pvalues(
            os.path.join(self.figuredir, "plot_moduleModel_") + self.data.paradict['output']
            )

    def draw_top_modules(self):
        '''
        May add later
        '''
        pass

    def export_top_modules(self):
        '''
        export to figures/network_modules/
        for visualizatin using CytoScape
        '''
        counter = 0
        for M in self.MA.top_modules:
            counter += 1
            M.export_network_txt(self.model, 
                                 os.path.join(self.moduledir, 'module_' + str(counter) + '.txt'))

    def export_activity_network(self):
            self.AN.export_network_txt(self.model, 
                                os.path.join(self.moduledir, 'activity_network.txt'))

    
    def export_cpd_attributes(self):
        '''
        Only for significant Trios.
        This exports attributes for network visualization in Cytoscape 3.
        
        mixedNetwork.hit_Trios is updated after initiation of ActivityNetwork
        [(mzFeature, EmpiricalCompound, cpd),...]
        
        '''
        s = 'CPD_ID\tCPD_NAME\tmz_row_Statistic\tFrom_mz_row\tFrom_EmpiricalCompound\n'
        for T in self.mixedNetwork.hit_Trios:
            name = self.model.dict_cpds_def.get(T[2], '')
            mz_row_Statistic = self.mixedNetwork.rowDict[T[0]].statistic
            s += '\t'.join([T[2], name, str(mz_row_Statistic), T[0], T[1].EID]) + '\n'
            
        with open(os.path.join(self.moduledir, "Node_attributes.txt"), 'w') as O:
            O.write(s)
    
    
    def export_network(self, network):
        '''
        Write ActivityNetwork for Cytoscape 3 formats.
        Changed to .txt in mummichog v 1.0.
        '''
        s = 'SOURCE\tTARGET\tENZYMES\n'
        for e in network.edges():
            s += e[0] + '\t' + e[1] + '\t' + self.network.MetabolicModel.edge2enzyme.get(e, '') + '\n'
        
        out = open(os.path.join(self.sifdir, 'ActivityNetwork.txt'), 'w')
        out.write(s)
        out.close()



class HtmlExport:
    '''
    This class is not tied to any data. 
    Serving as a skeleton to be used for export functions in AnalysisCentral.
    In future versions, this can use existing tools to convert JSON to HTML.
    '''
    def __init__(self):
        self.elements = []
        self.jsdata = ''
        
        self.HTML_HEAD = HTML_HEAD
        self.HTML_END = HTML_END
        self.javascript_HEAD = javascript_HEAD
        self.javascript_END = javascript_END
        
        
    def write_tag(self, s, tag, classname='', htmlid=''):
        cn = ''
        if htmlid: cn += ' id="' + htmlid + '"'
        if classname: cn += ' class="' + classname + '"'
        
        return '<' + tag + cn + '>' + s + '</' + tag + '>'


    def add_element(self, s, tag, classname='', htmlid=''):
        self.elements.append( self.write_tag(s, tag, classname, htmlid) + '\n')

    def make_select_menu(self, N):
        '''
        Showing
        <select id="nmodule_select" onchange="updateDrawing();">
          <option value="an">Activity network</option>
          <option value="module1">Module 1</option>
          <option value="module2">Module 2</option>
        </select>
        Node size
        <select id="nsize_select"  onchange="updateDrawing();">
          <option value=16>Medium</option>
          <option value=8>Small</option>
          <option value=32>Large</option>
        </select>
        '''
        s = '''Showing
        <select id="nmodule_select" onchange="updateDrawing();">
        '''
        for m in ["Module " + str(ii) for ii in range(1, N)]:
            s += "<option>%s</option>" %m
            
        return s + '''  <option>Activity network</option>
                   </select>
                    Node size <select id="nsize_select"  onchange="updateDrawing();">
                      <option value=16>Medium</option>
                      <option value=8>Small</option>
                      <option value=32>Large</option>
                    </select>
                    Visual style <select id="vstyle_select"  onchange="updateDrawing();">
                      <option value=0>Drag</option>
                      <option value=1>Force</option>
                    </select>
                    '''


    def make_js_data(self, networks, cpdnamedict, dict_cpd_statistic):
        '''
        Make data for both d3.js and cytoscape.js visualization.
        
        d3js visualization from multiple networks.
            var nodes = [
                [ {name:"Myriel",group:1},
                {name:"Napoleon",group:1},
                {name:"Mlle.Baptistine",group:1},
                {name:"Mme.Magloire",group:1},
                {name:"CountessdeLo",group:1},], ...];
                
                var links = [
               [ {source:1,target:0},
                {source:2,target:0},
                {source:3,target:0,value:10},
                {source:5,target:0,value:1},], ...];
        
        cytoscape.js visualization:
        	var cytonodes = [ [
				{ data: { id: "C01227", weight: 2.11 } },
				{ data: { id: "C05475", weight: 3.65 } },
				{ data: { id: "C05487", weight: 3.65 } },
				{ data: { id: "C03205", weight: -1.8 } },
					  ], ... ]
	  
					var cytoedges = [ [
				{ data: { id: "C01227C05498", weight: 1, source: "C01227", target: "C05498" } },
				{ data: { id: "C01227C05487", weight: 1, source: "C01227", target: "C05487" } },
				{ data: { id: "C01227C03205", weight: 1, source: "C01227", target: "C03205" } }, ... ]]
        
        Color coded by group. Thru additional argument colordict.
        The visualization can get a lot more complicated. Using only names for cytoscape.js for now. 
        When implementing other features, both id and name will be needed.
        
        '''

        cpdcolordict = self.rescale_color(dict_cpd_statistic)
        
        nodestr, edgestr, cynodestr, cyedgestr = '', '', '', ''
        for network in networks:
            nodestr += '['
            cynodestr += '['
            ordered_nodes = network.nodes()
            nodesdict, ii = {}, -1
            for n in ordered_nodes:
                ii += 1
                nodesdict[n] = ii
                nodestr += '{name:%s, group:%d}, ' %(quote(cpdnamedict.get(n, n).split(';')[0]), cpdcolordict.get(n, 10))
                cynodestr += '{ data: ' + '{id:%s, group:%d} }, ' %(quote(cpdnamedict.get(n, n).split(';')[0]), cpdcolordict.get(n, 10))
                
            edgestr += '['
            cyedgestr += '['
            for e in network.edges():
                edgestr += '{source:%d, target:%d}, ' %(nodesdict[e[0]], nodesdict[e[1]])
                cyedgestr += '{ data: { id: "%s", weight: 1, source: %s, target: %s } }, ' %('-'.join(e), 
                				quote(cpdnamedict.get(e[0], e[0]).split(';')[0]), quote(cpdnamedict.get(e[1], e[1]).split(';')[0]) )

            nodestr += '], '
            cynodestr += '], '
            edgestr += '], '
            cyedgestr += '], '
        
        total_d3_data = 'var nodes = [ ' + nodestr + '];\n\n        var links = [' + edgestr + '];\n\n'
        total_cytoscapejs_data = '        var cytonodes = [ ' + cynodestr + '];\n\n        var cytoedges = [' + cyedgestr + '];\n\n'
        self.jsdata = total_d3_data + total_cytoscapejs_data
        
    def rescale_color(self, dict_cpd_foldchange):
        '''
        rescale color to 0~20
        '''
        max_fc = max([abs(x) for x in dict_cpd_foldchange.values()])
        newdict = {}
        for k,v in dict_cpd_foldchange.items(): newdict[k] = int(10*v/max_fc) + 9
        return newdict
    
    def rescale_color_testing(self, dict_cpd_statistic):
        '''
        rescale color to 0~20
        7 steps, centered at 0, +- 3 * stdev
        '''
        stdev = np.std(dict_cpd_statistic.values())
        step = stdev
        newdict = {}
        for k,v in dict_cpd_statistic.items(): 
            if v > 3*step:
                newdict[k] = 9+9
            elif v < -3*step:
                newdict[k] = 1
            else:
                newdict[k] = int(v/step) + 9
        
        return newdict

    def export_text(self):
        s = ''
        for element in self.elements:
            s += element
            
        return [self.HTML_HEAD, 
                s + self.javascript_HEAD + self.jsdata + self.javascript_END, 
                self.HTML_END]



# --------------------------------------------------------
#
# a few function for utility or visualization scheme
#

def write_xlsx(data, filename, sheetname=''):
    '''
    data in regular lists, [ [], [], ... ]
    '''
    workbook = xlsxwriter.Workbook(filename, {'strings_to_numbers': True})
    ws = workbook.add_worksheet(sheetname)
    row = 0
    for line in data:
        col = 0
        for item in line:
            ws.write(row, col, item)
            col += 1
        row += 1
    
    workbook.close()

def write_yes_no_MassFeature(f):
    if f.is_significant:
        return "yes"
    else:
        return "no"

def quote(s):
    return '"'+s+'"'
