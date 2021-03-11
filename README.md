# Perturbome: Pipelines based on machine learning to identify core genes related to multiple perturbations
Pipelines for the analysis of the perturbome or response to multiple perturbations in Pseudomonas aeruginosa

Work: A first perturbome of Pseudomonas aeruginosa: Identification of core genes related to multiple perturbations by a machine learning approach

Developer: Jose A. Molina-Mora, PhD.

To run the three pipeline, just download them in a directory. After installing the required packages, you can run all the pipelines directly.
Comments are provided in each step.

Details of the files:
  1. Pipeline 1 Data & normalization: You can download data from GEO database and normalize using RMA (for expression values by Microarrays).
  2. Pipeline 2 Single Partition Method: a single partitioning strategy is done, and genes are ranked. Then performance of machine learning classifiers are tested and a topK values are selected. See Methods for details. 
  3. Pipeline 3 Multiple Partitions Method: multiple partitioning strategies are run with the help pf a special script "Classif_Independent.R", which is also provided. Consensus TopK genes are selected after multiple partitions are assessed. 
  4. Classif_Independent.R: Function necessary for Pipeline 3, in which the partitions, machile-learning classifiers and gene selection are executated. 
  5. Complete_Dataset: File with the normalized data (obtained from Pipeline 1), but manual selection of conditions (see Methods) was required to define the final list, as well as metadata was incorporates. Here you can find the GEO Ids by experiment and serie, name of the perturbation and other data.
  
Contact: 

Jose Arturo Molina Mora, Ph.D.
Microbiólogo bioinformático 
Tel. (+506) 2511 8646
e-mail: jose.molinamora@ucr.ac.cr
https://orcid.org/0000-0001-9764-4192

Sección de Biología Celular y Molecular, 
Facultad de Microbiología, 
Universidad de Costa Rica

