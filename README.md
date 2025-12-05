Mummichog 3 dev
===============

Mummichog is a Python program for analyzing data from high throughput, untargeted metabolomics.
It leverages the organization of metabolic networks to predict functional activity directly from feature tables,
bypassing metabolite identification. 

The last of version 2 is under branch mummichog-2.7.

This is version 3 under development.

Project is moved to new organization https://github.com/metabolomics-cloud, to follow examples of https://scverse.org/

## Input and Output

Input:
1. User supplied features with m/z, rtime, p-value from a statistical test. If no unique feature ID is supplied, row_number will be used as ID.
2. [Optional] Annotation table on the features. 
3. [Optional] Metabolic model to use. Default and optional models are provided by mummichog. Default in JSON (option for web app preload). 

Annotation can be from authentic standards and MS/MS. We use metDataModel to structure annotation. One can use JMS to perform annotation on a dataset. 

The format of a metabolic model contains list_of_pathways, list_of_reactions and list_of_compounds. 
We need chemical formulas of the compounds, which may not be available in a GSMM. 
- If compound formulas are provided in a model, be sure to have correct neutral mass and not a salt format. Salts are formed in both biological systems and in mass spectrometry. We use neutral formula to calculate adducts in mass spec, which includes salts. 
- If formulas are not provided in a model, we need to look them up via compound identifiers.
- The compound identifiers need to align with other data. This should be taken care of in mummichog supplied models. For developers, a translation module is needed. 

Some indexing and calculation on chemical formulas are involved in testing the match to metabolomic data patterns. Pathway definition may not be available in some models. "Subsystem" could be a substitute. 
Example of model conversion in:
- https://github.com/shuzhao-li/mummichog/blob/master/mummichog/models.py
- https://github.com/shuzhao-li-lab/JMS/tree/main/notebooks

Outpout are 
1. Result tables and figures, result.html as ver 2.
2. JSON strings from pathway analysis and network module analysis for programmatic use.

There's a separate repository for web-based mummichog tool, which handles UI and result visualization. 

## Planning

1. Moved 

2. New test datasets 

3. Milestone 3.1: Support of N metabolic models in JSON

4. Milestone 3.2: Run with new annotation formats, backward compatible and user-supplied annotation is optional

Azimuth DB should be renamed Mummichog DB.

Francisco and YC deployed web mummichog apps. Let's keep this as core package, with minimal dependency. 

## Test run

```
(base) MLG-JGM467:mummichog lish$ python3 -m mummichog.main -i tests/ineurons_ttest_1127.tsv -j testneuron -a tests/empCpds_with_annotations.json -d .
```

---
Old text -

## The mummichog suite includes

* mummichog(3): core algorithm package for pathway/network analysis

* cloud-mummichog: server and worker (RESTful) implementations

* Azimuth DB: the chemical database for biology, including metabolic models

* metDataModel: data models for metabolomics, used by mummichog and Azimuth DB

* mass2chem: common utilities in interpreting mass spectrometry data, annotation

* massBrowser: visualization using js (code reusable from CSM and mummichog web app)







## set up env for development (Python3, using virtualenv on Linux)

sudo apt install python3-dev python3-venv

sudo pip3 install virtualenv

virtualenv env

source env/bin/activate

A few libraries used for mummichog, for example:

(env) $ pip install scipy matplotlib xlsxwriter networkx

(env) $ deactivate

*To run test:*

shuzhao@canyon:~/li.github/mummichog3$ python3 -m mummichog.main -f mummichog/tests/testdata0710.txt -o t3


