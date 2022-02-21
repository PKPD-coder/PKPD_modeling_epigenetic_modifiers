# PKPD_modeling_epigenetic_modifiers

Read-me

This repository serves as an addendum to the following peer-reviewed article published in Pharmacology and Therapeutics (2022):

---------------------------------------------------------------------------------------------------------------------------
Leveraging Modeling and Simulation to Optimize the Therapeutic Window for Epigenetic Modifier Drugs

Antje Walz, Arthur J. Van De Vyver, Li Yu, Marc Birtwistle, Nevan J. Krogan, and Mehdi Bouhaddou

---------------------------------------------------------------------------------------------------------------------------

___Please cite accordingly.___

___BACKGROUND___

Despite a promising potential as an anti-cancer treatment, the clinical development of epigenetic modifier drugs is being hampered 
by the frequent occurence of dose-limiting toxicities of which thrombocytopenia is an important one. 
In the peer-reviewed article published in Pharmacology & Therapeutics, the authors discuss how modeling and simulation can be employed 
to study the pharmacology of epigenetic modifier drugs and how it can help to select optimal doses and dosing regimens to achieve optimal 
tumor growth inhibition while limiting the occurence of thrombocytopenia.

For the purpose of this article, a mathematical model was developed to describe the pharmacokinetics (PK) and 
pharmacodynamics (PD) of a hypothetical epigenetic modifier drug with a special focus on assessing the effects on tumor growth inhibition 
and thrombocytopenia by various dose levels and dosing regimens. 

This repository contains a Python file that allows the reader to replicate the simulations presented in the main article. 
The file contains a system of coupled Ordinary Differential Equations (ODEs) that describe the PK/PD profile of the drug. 
The reader can modify parameter values to test out hypothetical scenarios and explore ideal dosing regimens.

___MODEL___

The PK (described with a 2-compartmental model), biomarker response, and tumor growth inhibition (TGI) components of the model as 
well as the related parameter values were derived from the model developed by Bouhaddou et al. (1).The model component of drug 
effect on platelet formation is derived from the work from Chalret du Rieu et al. (2).

___THROMBOCYTOPENIA___

Suppression of circulating platelet counts that lead to varying grades of thrombocytopenia were assessed following the 
Common Terminology Criteria for Adverse Events (CTCAE) v5.0 hematologic criteria (3):

Grade 0 (baseline): 			150-450x10^9 cells/L (initial value used for simulation was 300x10^9 cells/L)
Grade 1 (low toxicity): 		75-150x10^9 cells/L
Grade 2 (moderate toxicity): 		50-75x10^9 cells/L
Grade 3 (severe toxicity): 		25-50x10^9 cells/L
Grade 4 (life-threatening toxicity): 	0-25x10^9 cells/L

When thrombocytopenia occurs over the course of treatment with an epigenetic modifier drug, the treating physician may decide to 
stop treatment until platelet counts sufficiently recover. Depending on the chosen dosing regimen, insufficient platelet recovery 
may occur before the next dose is given and these regimens may be suboptimal. We used a cut-off of 100x10^9 cells/L to determine 
whether re-dosing is allowed based on Piette and Broussard-Steinberg (4).

These considerations are taken into account in figure 4 and 5 of the manuscript, for which a multitude of different doses and 
dosing regimens were simulated. From each simulated dosing scenario, the achieved tumor growth inhibition as well as 
area-under-the-curve (AUC) of platelet counts over time were recorded. The prior allows us to rank the dosing scenarios 
in terms of efficacy, the latter to rank them in terms of toxicity.

Based on these two metrics, a figure of clinical utility can be generated (figure 4) where the AUC of platelets relative 
to baseline is given as the y-axis and tumor growth inhibition is given as the x-axis. 
Figure 4 enables us to differentiate the different dosing scenarios with respect to their safety and efficacy profiles.

Figure 5 is a companion figure to the clinical utility shown in figure 4. It shows how the platelet count AUC relative to baseline 
relates to the grade of thrombocytopenia achieved in the simulation (figure 5B) and whether there has been sufficient platelet recovery 
before the next dose was given (figure 5A).


___REFERENCES___

(1) Bouhaddou M, Yu LJ, Lunardi S, Stamatelos SK, Mack F, Gallo JM, et al. Predicting In Vivo Efficacy from In Vitro Data: 
Quantitative Systems Pharmacology Modeling for an Epigenetic Modifier Drug in Cancer. Clinical and translational science. 2020;13(2):419-29

(2) Chalret du Rieu Q, Fouliard S, White-Koning M, Kloos I, Chatelut E, Chenel M. Pharmacokinetic/Pharmacodynamic modeling of abexinostat-induced 
thrombocytopenia across different patient populations: application for the determination of the maximum tolerated doses in both lymphoma and solid 
tumour patients. Investigational new drugs. 2014;32(5):985-94.

(3) https://ctep.cancer.gov/protocoldevelopment/electronic_applications/docs/CTCAE_v5_Quick_Reference_8.5x11.pdf (accessed January 5th 2022)

(4) Piette WW, Broussard-Steinberg CM. 63 - Hematologic Toxicity of Drug Therapy. In: Wolverton SE, editor. 
Comprehensive Dermatologic Drug Therapy (Fourth Edition): Elsevier; 2021. p. 689-99.e4.
