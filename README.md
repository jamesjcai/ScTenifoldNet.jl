# scTenifoldNet.jl
 scTenifoldNet: Construct and Compare scGRN from Single-Cell Transcriptomic Data

This package provides an implementation of scTenifoldNet.
See [bioRxiv - scTenifoldNet: a machine learning workflow for constructing and comparing transcriptome-wide gene regulatory networks from single-cell data](https://doi.org/10.1101/2020.02.12.931469)
for more information.

Requires Julia 1.3 or higher.

## Installation

```
] add https://github.com/jamesjcai/scTenifoldNet.jl
```

## Usage

The simplest way is to just call the `transform` function with an array of numbers.

```
julia> using scTenifoldNet

julia> X0=rand(100,1000);
julia> X1=copy(X0)
julia> X1[4,:].=0.0
julia> Z0=tenrnet(X0)
julia> Z1=tenrnet(X1)
julia> d,aln0,aln1=manialn(Z0,Z1)
julia> fc,p,adjp=drgenes(d)

julia> using StatsPlots, Distributions
julia> x=rand(Chisq(1), length(fc)) 
julia> qqplot(x, fc)
```
Available functions:
--------------------

|Code| Function |
|:-|:-|
|pcnet|Computes a gene regulatory network based on principal component regression|
|tensordecomp|Performs CANDECOMP/PARAFAC (CP) Tensor Decomposition|
|manialn|Performs non-linear manifold alignment of two gene regulatory networks|
|drgenes|Evaluates gene differential regulation based on manifold alignment distances|
|tenrnet|Subsmple cells, construct single-cell gene regulatory networks (scGRNs) using principal component regression (pcnet), and denoise scGRNs using tensor decomposition (tensordecomp).|

Example:
--------
#### Loading scTenifoldNet
Once installed, **scTenifoldNet** can be loaded typing:
```{julia}
using scTenifoldNet
```

#### Simulating of a dataset 
Here we simulate a dataset of 2000 cells (columns) and 100 genes (rows) following the negative binomial distribution with high sparsity (~67%).
```{julia}
d=NegativeBinomial(20,0.98)
X=rand(d,100,2000)
```

#### Generating a perturbed network 
We generate a perturbed network modifying the expression of genes 10, 2, and 3 and replacing them with the expression of genes 50, 11, and 5.
```{julia}
Y=copy(X)
Y[10,:]=Y[50,:]
Y[2,:]=Y[21,:]
Y[3,:]=Y[5,:]

X=X[:,vec(sum(X,dims=1).>30)]
Y=Y[:,vec(sum(Y,dims=1).>30)]
```
#### scTenifoldNet
Here we run **scTenifoldNet** under the H0 (there is no change in the regulation of the gene) using the same matrix as input and under the HA (there is a change in the regulation of the genes) using the control and the perturbed network.
```{julia}
Z0=scTenifoldNet.tenrnet(X, donorm=true)
Z1=scTenifoldNet.tenrnet(Y, donorm=true)
```
#### Differential regulation based on manifold alignment distances
As is shown below, under the H0, none of the genes shown a significative difference in regulatory profiles using an FDR cut-off of 0.1, but under the HA, the 6 genes involved in the perturbation (50, 11, 2, 10, 5, and 3) are identified as perturbed.
```{julia}
d,aln0,aln1=scTenifoldNet.manialn(Z0,Z1)
fc,p,adjp=scTenifoldNet.drgenes(d)
```

#### Plotting the results
Results can be easily displayed using quantile-quantile plots. Here we labeled in red the identified perturbed genes with FDR < 0.05.
![Example](https://raw.githubusercontent.com/cailab-tamu/scTenifoldNet/master/inst/readmeExample.png)
```{julia}
using StatsPlots
x=rand(Chisq(1), length(fc)) 
qqplot(x, fc)
```

Citation
--------
To cite **scTenifoldNet** in publications use:

  Daniel Osorio, Yan Zhong, Guanxun Li, Jianhua Huang and James Cai (2019). scTenifoldNet: Construct and Compare scGRN from Single-Cell Transcriptomic Data. R package version 1.2.0.
  https://CRAN.R-project.org/package=scTenifoldNet

A BibTeX entry for LaTeX users is
```
  @Manual{,
    title = {scTenifoldNet: Construct and Compare scGRN from Single-Cell Transcriptomic Data},
    author = {Daniel Osorio and Yan Zhong and Guanxun Li and Jianhua Huang and James Cai},
    year = {2019},
    note = {R package version 1.2.0},
    url = {https://CRAN.R-project.org/package=scTenifoldNet},
  }
  ```

