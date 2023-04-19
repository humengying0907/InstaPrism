# InstaPrism
**InstaPrism** is R package to deconvolute cellular proportion and gene expression in bulk RNA-Seq data. 
It provides fast deconvolution implementation of [BayesPrism](https://github.com/Danko-Lab/BayesPrism) while maintaining highly comparable performance. 
## Installation
```````
library("devtools");
install_github("humengying0907/InstaPrism")
```````
## Running time comparsion with BayesPrism
Below is a running time comparsion when running deconvolution on the tutorial [data](https://github.com/Danko-Lab/BayesPrism/tree/main/tutorial.dat) 
provided in BayesPrism. 

![fig4](https://user-images.githubusercontent.com/54827603/219739883-1b5c4f3b-e2cb-4843-a9d1-10e0381b17f8.png)



## Tutorial
Check [InstaPrism_tumorial](https://humengying0907.github.io/InstaPrism_tutorial.html) for detailed implementation of InstaPrism and compare its performance with BayesPrism.

## Reference
M. Hu and M. Chikina, “InstaPrism: an R package for fast implementation of BayesPrism.” bioRxiv, p. 2023.03.07.531579, Mar. 12, 2023.
doi: https://doi.org/10.1101/2023.03.07.531579
