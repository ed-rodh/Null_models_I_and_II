# Null_models_I_and_II

 Simulating annealing based strategies to randomize sequences data composed by:
 
  + Null model I: ramdomizes the input MSA by conserving the single columns statistics.
 
  + Null model II: randomizes the input MSA preserving both single columns statistics and pairwise hamming distances between sequences.
   
 A complete description of the algorithm can be found at (Edwin Rodriguez Horta, Weigt M.,Phylogenetic correlations have limited effect on coevolution-based contact prediction in proteins).

The code now requires at least Julia version 1.0 or later.

# Install
To install the package under version >= 1.0 use the package manager 

(v1.?) pkg> add https://github.com/ed-rodh/Null_models_I_and_II

# Overview


# Usage
Load the julia package by

 ```julia
 julia> using Null_models_I_and_II
 ```
This software provide two main functions ```sample_from_Null_model_I(infile::String;outfile="",Temp=10.0)``` and ```sample_from_Null_model_II(infile::String;outfile="",schuffle=1,Temp=10.0)``` which are in charge of generating new sequence data according to null models I and II respectively. Function ```julia sample_from_Null_model_I()``` is faster. 

Both functions take as input the name of multiple sequence alignments in fasta format. Optional parameter   "outfile"    name text file for the output alignment exported in fasta format and parameter "Temp" define the  temperature for  shuffling alignment via Markov chain Monte Carlo (MCMC).   


# Output

# Requirements

 
