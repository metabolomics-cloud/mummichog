Mummichog
=========

Mummichog is a Python program for analyzing data from high throughput, untargeted metabolomics.
It leverages the organization of metabolic networks to predict functional activity directly from feature tables,
bypassing metabolite identification. The features include

* computing significantly enriched metabolic pathways
* identifying significant modules in the metabolic network
* visualization of top networks in web browser
* visualization that also plugs into Cytoscape
* tentative annotations
* metabolic models for different species through plugins

This is code branch for mummichog package version 2.7. 
Master branch is on version 3.

Mummichog version 2.2 was the last version using Python 2; in branch as mummichog-python2.

Installation
------------

Mummichog can be installed using pip (pip Installs Packages), the Python package manager. The command below will install the default (version 2):

    pip install mummichog

This is OS independent. To read more on [pip here](https://pip.pypa.io/en/stable/installing/#installing-with-get-pip-py).

One can also run mummichog without installing it. Direct python call on a downloaded copy can work, e.g.

    python3 -m mummichog.main -f mummichog/tests/testdata0710.txt -o t2


Use custom metabolic model
--------------------------

The `-n` argument now (v2.6) takes user specified metabolic model in JSON format, e.g.

    python3 -m mummichog.main -n mummichog/tests/metabolicModel_RECON3D_20210510.json -f mummichog/tests/testdata0710.txt -o t3

The porting of metabolic model is demonstrated 
[here https://github.com/shuzhao-li/Azimuth/blob/master/docs/](https://github.com/shuzhao-li/Azimuth/blob/master/docs/From-GEM-to-metDataModel-20210510.ipynb)

Please note that identifier conversion is a major issue in genome scale models. Users benefit greatly from including chemical formula (neutral_formula) and molecular weight (neutral_mono_mass) in the model.

Citation
--------
The initial paper on mummichog is described in Li et al. Predicting Network Activity from High Throughput Metabolomics. PLoS Computational Biology (2013); doi:10.1371/journal.pcbi.1003123. 

