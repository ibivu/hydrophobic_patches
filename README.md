## Abstract

Motivation: Proteins tend to bury hydrophobic residues inside their core during the folding process to provide stability to the protein structure and to prevent aggregation. Nevertheless, many proteins do expose such 'sticky' hydrophobic residues to the solvent. Hydrophobic residues may play an important functional role, for example in protein-protein interactions, ligand binding and interactions with the membrane. Here, we investigate how hydrophobic (sticky) the surface of a protein is by providing measures for surface hydrophobicity, and by predicting these measures from sequence. Finally, we analyse the over-expression of sticky proteins within the human proteome by considering tissues, pathways and diseases in which such proteins occur.
Results: Firstly, we define three structure-based measures: the total hydrophobic surface area (THSA), the relative hydrophobic surface area (RHSA) and - using our MolPatch method - the largest hydrophobic patch (LHP). Secondly, we analyse how easy it is to predict these measures from sequence: by adapting solvent accessibility predictions from NetsurfP2.0, we obtain well-performing prediction methods for the THSA and RHSA, while the LHP is more difficult to predict. We show that very hydrophobic proteins are typically not highly expressed, suggesting there is evolutionary pressure against overabundant sticky proteins. Despite this, we show that sticky proteins with large hydrophobic surface areas are surprisingly abundant in the human brain, kidneys and blood cells and that such proteins are associated with neurodegenerative pathways.

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
