# Laboratory Bioinformatic-1 Project
## written by: Behrouz Mollashahi, Department of Pharmacy and Biotechnology (FaBiT), University of Bologna, 40126 Bologna, Italy
### This project is written for the laboratoy bioinformatic 1 course, Bologna university
## Abstract
Introduction: Kunitz domain is a serin protease domain exerting various functions within the cells. They are widely studied in plants, animals, and humans. Various studies show their importance in various fields. The aim of this study is to build a Hidden Markov Model (HMM) that can
predict the Kunitz domain and evaluate the model.
Method: In this study, data was extracted from Protein Data Bank (PDB) website under constraints and then Multiple structural Alignment (MSA) was performed with the PDB-fold website, the results were then cleaned and used to generate HMM model with the HMMER package. Finally, the evaluation of the model was done with the accuracy test and Matthew Correlation Coefficient (MCC).
Result: Based on testing the generated HMM model with the Swiss-prot database the high MCC score shows a significant performance of the HMM model generated for the Kunitz domain.
