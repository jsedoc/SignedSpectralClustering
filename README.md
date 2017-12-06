# SSC

A good place to start to understand the code is from example_clusters.m. I still
need to add further documentation. ncutsK_v5 is interactive code. This is largely
for visualization and inspection.

The main function for usage is sncut(W,K) where W is the weight matrix and K
is the number of clusters. 

Optional type ... sncut(W,K,'type','hard'||'flex'||'soft')

Please send questions to Joao Sedoc (joao at cis dot upenn dot edu).

For the theory behind SSC see :

Jean Gallier ,"Spectral Theory of Unsigned and Signed Graphs. Applications to Graph Clustering: a Survey" https://arxiv.org/abs/1601.04692

Sedoc et al., "Semantic Word Clusters Using Signed Normalized Graph Cuts" https://arxiv.org/abs/1601.05403


# TBD: I will move the code for the paper to be under the folder synonym_clusters

