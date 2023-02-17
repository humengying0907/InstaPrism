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

![fig4](https://user-images.githubusercontent.com/54827603/219703956-f56c88e0-e1ca-49ef-b76d-115eb9ec1da5.png)


## Tutorial
Check the [InstaPrism_tumorial](https://humengying0907.github.io/InstaPrism_tutorial.html) for details.
