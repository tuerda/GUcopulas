# GUcopulas
Code to simulate from posteriors of grid-uniform copula models in julia


Example code:
```
include("mcmccopula.jl")
using Distributions

#data. Each data point should be a d dimensional vector. from a copula

function makedat() #This creates data from a (half & half) mixture of gaussian copulas.
    normies=rand(MvNormal([0,0],[1 0.5; 0.5 1]),1000)
    normies2=rand(MvNormal([0,0],[1 -0.9; -0.9 1]),1000)
    ddat=hcat(copy(normies),copy(normies2))
    for (i,norm) in enumerate(eachcol(ddat))
        ddat[:,i]=[cdf(Normal(0,1),norm[1]), cdf(Normal(0,1),norm[2])]
    end
    return ddat
end
dat=makedat()

##Lets do it!!

startfixed=fixedState(2)
inferfixed=Fixedchain(startfixed,[],dat)
mcmc(inferfixed,15000)

#Copulas stored in inferfixed.history

```
