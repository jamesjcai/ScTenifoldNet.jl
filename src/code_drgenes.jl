# ] add https://github.com/tk3369/BoxCoxTrans.jl
# ] add UnicodePlots
using Distributions, UnicodePlots, BoxCoxTrans, MultipleTesting

"""
x = rand(Gamma(2,2), 10000) .+ 1;
histogram(x)
histogram(BoxCoxTrans.transform(x))
BoxCoxTrans.lambda(x).value
histogram(BoxCoxTrans.transform(x, 0.01))
histogram(BoxCoxTrans.transform(x; scaled = true))

scatterplot(randn(50), randn(50), title = "My Scatterplot")
plt = lineplot([cos, sin], -π/2, 2π)
lineplot!(plt, -0.5, .2, name = "line")
barplot(["Paris", "New York", "Moskau", "Madrid"],
        [2.244, 8.406, 11.92, 3.165],
        title = "Population")
"""

using Statistics, StatsBase
d=rand(30)
zscore(d)
(d.-mean(d))./std(d) == zscore(d)
d²=d.^2
FC=d²./mean(d²)

1-cdf(Chisq(1),7.8)==ccdf(Chisq(1),7.8)



pvals = ccdf.(Chisq(1),FC)
pAdjusted = MultipleTesting.adjust(pvals, BenjaminiHochberg())


# https://stats.stackexchange.com/questions/171074/chi-square-test-why-is-the-chi-squared-test-a-one-tailed-test

1-cdf(Chisq(1),4) == 2*(1-cdf(Normal(0,1),2))
ccdf(Chisq(1),4) ≈ 2*(1-cdf(Normal(0,1),2))

#pValues <- pchisq(q = FC,df = 1,lower.tail = FALSE)
#pAdjusted <- p.adjust(pValues, method = 'fdr')

using Distributions
# https://rosettacode.org/wiki/Verify_distribution_uniformity/Chi-squared_test#Julia 
function eqdist(data::Vector{T}, α::Float64=0.05)::Bool where T <: Real
    if ! (0 ≤ α ≤ 1); error("α must be in [0, 1]") end
    exp = mean(data)
    chisqval = sum((x - exp) ^ 2 for x in data) / exp
    pval = ccdf(Chisq(2), chisqval)
    return pval > α
end
 
data1 = [199809, 200665, 199607, 200270, 199649]
data2 = [522573, 244456, 139979,  71531,  21461]
 
for data in (data1, data2)
    println("Data:\n$data")
    println("Hypothesis test: the original population is ", (eqdist(data) ? "" : "not "), "uniform.\n")
end
