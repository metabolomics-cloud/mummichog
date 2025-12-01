# About annotation for mummichog

Annotation can be from one or more methods, but often of varying confidence. 

This creates a complexity in modeling annotation. The ideal input to mummichog is a list of identities with probabilities.  

It's not a trivial problem. Examples are inclded in notebooks to illustrate how:
- khipu does pre-annotation in JSON
- khipus are chained to include MS/MS
- then to include LC-MS authentic library

Most features still have no annotation. 
Mummichog handles those by internal matching to metabolic models, considering annotation from other sources (provided by users). 
Using this approach, targeted metabolomics data become a special case of unambiguous annotation.

In mummichog, we try to field multiple situations where users supply no or parital annotations.
