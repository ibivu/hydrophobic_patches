# Abstract

In general, hydrophobic residues tend to be buried in the interior of a protein to avoid contact with a hydrophilic solvent. Nevertheless, large hydrophobic areas can be found on a protein’s surface. The size of the hydrophobic surface area and, more specifically, hydrophobic patches play essential roles in protein-protein interactions or non-specific aggregation. In the absence of structural information, the prediction of the size of the solvent-accessible hydrophobic area and patches is, therefore, valuable for deciphering its function. In this work, we present simple methods for predicting the total hydrophobic surface area and relative hydrophobic surface area of a protein. The methods were benchmarked against NetSurfP2. We show that solely NetSurfP2 results are a good estimator for the hydrophobic surface area. Furthermore, we present MolPatch: a method to calculate the hydrophobic patches of a protein from structural data. The three largest hydrophobic patches obtained with MolPatch contained a significant increase in protein interaction sites compared to random patches. We also show that the size of the largest patch can be predicted from the sequence with reasonable accuracy. These predictions can help with the identification of protein-protein binding sites and can lead to a better understanding of biological effects elicited by hydrophobic features.

## Project

This repository has five main directories.

* /data
  * all the generated csv files (will be created when running the code)
* /model
  * all the prediction models as pickles
* /research
  * jupyter notebook with figures for manuscript
* /src
  * the source code of this project

## Contact

dr. Sanne Abeln, s.abeln@vu.nl
MSc Juami van Gils, j.h.m.van.gils@vu.nl
