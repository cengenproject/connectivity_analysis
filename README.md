# Connectivity analysis through network differential gene expression (nDGE)
Hypothesis testing the genetic basis of the connectome. To run the code, modify the inputs in "nerve_ring_analysis.m" to select the desired neuron or neuron pair and press run.

A) Overview of standard differential gene expression analysis (DGE): DGE can be represented as a regression model with the response variables as the neuron-wise gene expression data and the design matrix as 1s and -1s denoting the group memberships of neurons. B) Schematic of network differential gene expression analysis (nDGE): nDGE is a generalization of DGE that represents the pairwise co-expression of all genes in all pairs of neurons as the response variable. The design matrix is a square matrix of all pairs of neurons with 1s and -1s placed in locations where pairs of edges of the network are contrasted. C) Design matrix for nDGE that yields identical results as standard DGE. D) Design matrix for nDGE that explores the global genetic differences between synapses and non-synaptic membrane contact. E) Design matrix for nDGE that explores the differential gene co-expression differences between the synapses of two different neurons. F) Design matrix for nDGE that aims to demonstrate gene co-expression enrichment in synapses of a particular neuron compared to non-synaptic membrane contacts. G) The statistical significance of nDGE test statistic is obtained by computing the mean and variance of the null distribution of the test statistic resulting from pseudoconnectomes that respect the network topology of the original connectome.


![overview](https://github.com/cengenproject/connectivity_analysis/blob/main/method.png)


