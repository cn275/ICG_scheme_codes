# ICG_scheme_codes
Fitting a coarse-grained model onto a fully atomistic model

Coarse-grained (CG) models are simplified representations of fully atomistic models as seen in MD, DFT, and MC. They often are developed to allow for much faster and larger simulations of a molecule of interest and have been used for simulating various properties (thermal, mechanical, phase stability, etc.) of interest. 

Due to the increased agility and ease of CG models some studies use them from the onset. Studies of this kind have shown unique/desirable behaviors/properties with CG models but lack the atomistic representation that will behave like these models. 

We have developed a methodology for taking an atomistic model and finding several mappings of the atomistic structure onto the coarse-grained structure of interest. Small simulations of the atomistic model can be used to calculate how well the atomistic model "fits" the CG model and be compared to how well other candidate molecules fit the CG model (https://doi.org/10.1021/acs.jcim.9b00232). 

This repository contains the code used in this study for generating the mappings and will be updated with updated codes. In particular a new version of the code will be developed using python and commonly used modules.
