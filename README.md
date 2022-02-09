# Inverse_Theory_Bayes_Codes

This repository is a compendium of codes I made during the course "Inverse Problem Theory - A Bayesian Approach".  This course was given by Prof. Klaus Mosegaard, during his visit at ETH Zurich in I-2017.

The most interesting results are in the solutions of assignements 4 and 5, where I solved the inversion problem given the data of an horizontal slice of 3D seismic reflection data.

I obtain statistics (histogram and marginal probabilities scanning all possible configurations) from a training image by using a 2x2 moving window. 

Then, i)an initial model is created and ii)conditional probabilities are computed based on the distribution (statistics) of the training image. This allows building a new image by picking random configurations from this distribution.

Finally, I solve the inversion problem using the Metropolis algorithm. I proposed a model which explain the observed data within the error bars.
