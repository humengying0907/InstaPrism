# InstaPrism
**InstaPrism** is R package to deconvolute cellular proportion and gene expression in bulk RNA-Seq data. 
It provides fast deconvolution implementation of [BayesPrism](https://github.com/Danko-Lab/BayesPrism) while maintaining same performance. 
## Installation
```````
library("devtools");
install_github("humengying0907/InstaPrism")
```````
## Running time comparsion with BayesPrism
Below is a running time comparsion when running deconvolution on the tutorial [data](https://github.com/Danko-Lab/BayesPrism/tree/main/tutorial.dat) 
provided in BayesPrism. 

![rt](https://user-images.githubusercontent.com/54827603/217444704-027b794f-1ac0-42c0-9bd0-e9a5d1bb338e.png)

## Tutorial
Check the [InstaPrism_tumorial](https://github.com/humengying0907/InstaPrism/blob/master/InstaPrism_tutorial.html) for details. \
To view, download the html file or paste the raw html to any [html viewer](https://codebeautify.org/htmlviewer).
