# InstaPrism
**InstaPrism** is R package to deconvolute cellular proportion and gene expression in bulk RNA-Seq data. Based on the same conceptual framework and corresponding generative mode from [BayesPrism](https://github.com/Danko-Lab/BayesPrism), **InstaPrism** re-implements BayesPrism in a derandomized framework by replacing the time-consuming Gibbs sampling steps in BayesPrism with a fixed-point algorithm, which greatly accelerated the calculation speed while maintaining highly comparable performance.
## Installation
```````
library("devtools");
install_github("humengying0907/InstaPrism")
```````
## Unique Features
Explore the following distinctive features of our package:

* Computational and Memory Efficiency: Optimized algorithms for efficient utilization of computational resources and memory
* Progress Bar Monitoring: Real-time progress monitoring of the deconvolution process
* Convergence Check: Summary table and plot to monitor algorithm convergence status
* Selective Reference Update: Ability to update the reference specifically for desired cell types of interest
* Cell-Type Specific Expression Reconstruction: Ability to reconstruct cell-type specific expression for bulk samples
## Running time comparsion with BayesPrism
Below is a running time comparsion when running deconvolution on the tutorial [data](https://github.com/Danko-Lab/BayesPrism/tree/main/tutorial.dat) 
provided in BayesPrism. 

<img src="https://github.com/humengying0907/InstaPrism/assets/54827603/7ecf1817-e04b-4739-994a-fabf3a1064d1" width=45% height=45%>


## Tutorial
Check [InstaPrism_tumorial](https://humengying0907.github.io/InstaPrism_tutorial.html) for detailed implementation of InstaPrism and compare its performance with BayesPrism.

## Reference
M. Hu and M. Chikina, “InstaPrism: an R package for fast implementation of BayesPrism.” bioRxiv, p. 2023.03.07.531579, Mar. 12, 2023.
doi: https://doi.org/10.1101/2023.03.07.531579
