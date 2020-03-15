cd(dirname(@__FILE__))

include("scTenifoldNet.jl")
using .scTenifoldNet
X0=rand(100,1000);
X1=copy(X0)
X1[4,:].=0.0

@show Threads.nthreads()

@time Z0=scTenifoldNet.tenrnet(X0, donorm=false)
@time Z1=scTenifoldNet.tenrnet(X1, donorm=false)
@time d,aln0,aln1=scTenifoldNet.manialn(Z0,Z1)
fc,p,adjp=scTenifoldNet.drgenes(d)

#using StatsPlots, Distributions
#x=rand(Chisq(1), length(fc)) 
#qqplot(x, fc)