# Conservation- and phylogeny-aware Null models for protein MSA

 Codes implement simulated-annealing based strategies to randomize sequences data, keeping conservation (Null model I) and phylogenetic (Null model II) information, but removing any signals of functional or structural coevolution. 
 
A complete description of the algorithm can be found in

Edwin Rodriguez Horta & Martin Weigt,Phylogenetic correlations have limited effect on coevolution-based contact prediction in proteins, xxx.

Please cite this paper in any publication using our codes.

The codes are written in Julia (https://julialang.org/ ) and require Version 1.0 or later.

# Install
To install the package under Julia Version >= 1.0 use the package manager (to be activated from the REPL using the key ] )

(v1.?) pkg> add https://github.com/ed-rodh/Null_models_I_and_II

# Overview
There are two strategies:
 
  + Null model I: ramdomizes the input multiple sequence alignments (MSA) by preserving the single columns statistics. The null model is thought to preserve conservation-related information, but to destroy any correlation due to phylogeny or residue coevolution.
  + Null model II: randomizes the input multiple sequence alignments (MSA) preserving both single columns statistics and pairwise hamming distances between sequences. The null model is thought to preserve conservation- and phylogeny-related information, but to destroy any correlation due to residue coevolution.
   
The code runs on MSA of amino-acid sequences in fast format. It assume that each entry of the MSA is a letter in { A, C, ..., Y,-} being either one of the 20 amino acids or the alignment gap “–”. The total alphabet size is q = 21.

# Usage
After installation, load the Julia package by

 ```julia
 julia> using Null_models_I_and_II
 ```
This code provide two main functions: ```sample_from_Null_model_I(infile::String;outfile="",shuffle_temp=10.0)``` and ```sample_from_Null_model_II(infile::String;outfile="",shuffle=1,shuffle_temp=10.0,T_factor_slow=0.8,T_factor_fast=0.1,min_temp="default",num_iter_max=20000000)``` which are in charge of generating new sequence data according to null models I and II as described above. Function ```sample_from_Null_model_I()``` is faster while ```sample_from_Null_model_II()``` may become very slow for deep MSA.. 

 + Both functions take as input a file containing a protein MSA in fasta format(parameter "infile"). 
 + The optional parameter "outfile" gives the name for the text file, where the output MSA is exported to.
 + The parameter "shuffle_temp" stands for the formal temperature in Monte Carlo sampling to shuffling the sequences data (is recommended to set shuffle_temp>1).
 + The parameter 'shuffle' allows (shuffle=1) or not (shuffle=0) to initialize the sampling strategy of Null model II from sequences data generated by Null model I. The use of shuffle=1 is strongly recommended to remove all functional and structural coevolutionary signals from the MSA before re-establishing the phylogenetic signal.
 
For function ```sample_from_Null_model_II()``` a gradual reduction of the temperature toward minimal temperature ('min_temp') is carried out by a linear annealing schedule with two slopes characterized by parameters 'T_factor_fast' and 'T_factor_slow'.
   + 'T_factor_fast' represent the rate of temperature reduction for high temperatures where the system is expected to wander initially towards a broad region of the search space (acceptance rate high).  
   + 'T_factor_slow'  represent the cooling rate at low temperatures where the search space become narrower (acceptance rate small).
   + This code assume a default 'min_temp' parameter  linked to the number of sequences in the MSA: min_temp=2/(10*M).
   + The parameter num_iter_max stand for the maximum number of swaps attempted.    

We also provide other potentially useful functions:
 + ```translate_fasta_to_num_matrix(msa_fasta::String,...)``` to read the MSA as fasta file and transform to an MSA as a matrix of numbers. The mapping is {-, A, C, ..., Y} -> {1, 2, 3,..., 21}.
 + ```transform_MSA_fasta(msa::Array{Int64,2},...)``` transforms an MSA from numbers to letters, using the inverse mapping {1, 2, 3,..., 21} -> {-, A, C, ..., Y}.
 + ```export_fasta_file(outfile::String,msa::Array{Char,2})``` export a MSA in fasta format from a matrix of letters.
 + ```PairwiseHammingDist(sample::Array{Int64,2})``` computes all pairwise Hamming distances between sequences from a numeric MSA stored in the matrix "sample".
# Output

The functions output a type Array{Char,2}, which describes a matrix of aminoacids (letters). If “outfile” isn’t empty the output MSA is exported to the indicated text file in fasta format.

=> For a quick test we recommend the alignment included in Null_models_I_and_II/data_test
# Requirements

The code requires at least Julia version 1.0 or later.

