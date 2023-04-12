##
## MCMC over grid-uniform copulas. Fixed grid version.
## Fixed grid works. Free grids in another file don't work well yet.
##

#Import copula objects
include("copulaug.jl") #Brings in Plots already.

##
##### parameters. In this case defining the parameter space requires a fair amount of code
##

function pixelate(G1vals, divs=10) #changed to add divs.
    pixels=zeros(0,divs)
    for i in 1:divs
        line=[]
        for j in 1:divs
            push!(line,sum(G1vals[((i-1)*Integer(floor(100/divs))+1):(i*Integer(floor(100/divs))),((j-1)*Integer(floor(100/divs))+1):(j*Integer(floor(100/divs)))]))
        end
        pixels=vcat(pixels,line')
    end
    return pixels
end

#Takes a bivariate copula density and a number of divisions and creates a 2 dimensional equally spaced grid uniform copula.
function makegucopfromfn(f,numdivs)
    grid=(1/numdivs):(1/numdivs):1
    divs=[grid,grid]
    sizes=((1/numdivs)^2)*ones(numdivs,numdivs)
    G1space=0.005:0.01:1
    G1vals=zeros(100,100)
    for (i,g0) in enumerate(G1space)
        for (j,g1) in enumerate(G1space)
            G1vals[i,j]=f([g0,g1])/10000
        end
    end
    probs=pixelate(G1vals, numdivs)
    return Copula(2,divs,sizes,probs)
end

#The distance function for the prior.
#Different parameter values are used to specify specific cases.
#As a result there are some discontinuities with respect to a and gamma.
function curlyD(cop::Copula, C0::Copula, a::Float64, gamma::Float64) #L2 Norm of difference to G0. If a=0 it just returns 1
    if a==0 #Prior is flat. Calculate nothing.
        return 1
    end
    if gamma==0 #No smoothing, use L2 norm.
        copdens=cop.probs ./ sqrt.(cop.sizes)
        copdensq=C0.probs ./ sqrt.(C0.sizes)
        return sum((copdens .- copdensq).^2)
    end
    R=CartesianIndices(cop.probs)
    Ifirst,Ilast=first(R),last(R)
    I1=oneunit(Ifirst)
    result=similar(cop.probs)
    cminusc0=cop.probs .- C0.probs
    if gamma==1 #ICAT case.
        for I in R
            NB=cminusc0[max(Ifirst, I-I1):min(Ilast, I+I1)]
            result[I]=sum((cminusc0[I].-NB).^2)
        end
        return sum(result)
    else #PCAT case.
        for I in R
            NB=cminusc0[max(Ifirst, I-I1):min(Ilast, I+I1)]
            absnb=length(NB) - 1 #a cell is not its own neighbor.
            sumand1=absnb * cminusc0[I]^2
            sumA=sum(NB)
            sumA -= cminusc0[I] #a cell is not its own neighbor
            result[I]= sumand1 - gamma * cminusc0[I] * sumA
        end
        return sum(result)
    end
end


##
## ##########The structure of a single step in the MCMC chain##########
##

abstract type AbstractState end

mutable struct fixedState <: AbstractState
    cop::Copula #The copula corresponding to the state.
    C0::Copula #Reference copula for the prior.
    alpha::Float64 #alpha parameter for the prior.
    energy::Float64 #energy of the copula
    pri::Float64 #value of the prior
    lik::Float64 #likelihood of the copula
    gamma::Float64 #gamma parameter for the prior
end

#Constructor Some variations are calculated automatically. If you have something more complex
#than what is done automatically, you can also initialize the object by hand.
function fixedState(dim=2; grid=Array(0.1:0.1:1), gamma=0., alpha=0., g0=x->1)
    divs=[]
    for a in 1:dim
        push!(divs,grid)
    end
    bits=ones(Int, dim) .* length(grid)
    sizes=((grid[1])^dim)*ones(bits...)
    probs=copy(sizes)
    cop=Copula(dim,divs,sizes,probs)
    C0=makegucopfromfn(g0, length(grid))
    prior=logprior(cop, alpha, gamma, C0)
    llik=0
    engy=prior+llik
    return fixedState(cop, C0, alpha, engy, prior, llik, gamma)
end


##########The necessary functions to get State off the ground##########

function loglik(data, cop::Copula; cache=[], remembered=true) #log likelihood function
    if data==[] #no data, simmulating from the prior.
        return 0
    end
    if cache==[]
        lik=0
        for d in eachcol(data)
            if ~ remembered
                here=zeros(Int,cop.dim)
                for (i,c) in enumerate(d)
                    locingrid=1
                    while (locingrid<length(cop.divs[i]))&(cop.divs[i][locingrid]<c)
                        locingrid+=1
                    end
                    here[i]=locingrid
                end
            else
                here=d
            end
            lik+=log(cop.probs[here...])-log(cop.sizes[here...])
        end
        return lik
    else
        return sum(log.(cop.probs ./ cop.sizes) .* cache)
    end
end

function logprior(cop::Copula, a=0., gamma=0., C0=cop) #Log-prior
    return (-(a/2)*curlyD(cop,C0, a, gamma))
end

# Kernels

function recexg(st::AbstractState, dats; cach=[], memorized=true) #make a move using the rectangle exchange kernel
    probsprop=rectangle(st.cop)[1]
    prop=deepcopy(st.cop)
    prop.probs=probsprop
    pri=logprior(prop, st.alpha, st.gamma, st.C0)
    lik=loglik(dats, prop, cache=cach, remembered=memorized)
    en=pri+lik
    logMH=min(0,en-st.energy)
    if log(rand())<logMH
        st.cop=prop #Accept
        st.energy=en
        st.pri=pri
        st.lik=lik
    end
end


##
## ##########MCMC time!##########
##

abstract type CopChain end

mutable struct Fixedchain <:CopChain #When the state is fixed
    now::fixedState
    history::Array{fixedState}
    data
end

mutable struct Cachedchain <:CopChain #When the state is fixed adding a cache to make likelihood calculation faster. This may or may not speed things up depending on the sample size, grid size, etc.
    now::fixedState
    history::Array{fixedState}
    data::Array{Float64}
    cache::Array{Float64}
end

function makecache!(ch)
    if ch.data==[]
        print("No data; this will simulate only from the prior.")
        return
    end
    divs=ch.now.cop.divs
    ch.cache=zeros(size(ch.now.cop.sizes)...)
    for d in eachcol(ch.data) # each datum
        here=zeros(Int, ch.now.cop.dim)
        for (i,c) in enumerate(d) # each dimension in the datum
            locingrid=1
            while (locingrid<length(divs[i])) & (ch.now.cop.divs[i][locingrid]<c)
                locingrid+=1
            end
            here[i]=locingrid
        end
        ch.cache[here...]+=1
    end
end

function memorize(data, cop::Copula)
    if all(a>=1 for a in data)
        return data
    end
    heres=zeros(Int, size(data))
    for (i,d) in enumerate(eachcol(data))
        here=zeros(Int,cop.dim)
        for (i,c) in enumerate(d)
            locingrid=1
            while (locingrid<length(cop.divs[i]))&(cop.divs[i][locingrid]<c)
                locingrid+=1
            end
            here[i]=locingrid
        end
        heres[:,i]=here
    end
    return heres
end


function mcmc(ch::Cachedchain, numberofsteps=15000; thinning=1, silent=false) #All of the magic happens here. Fixed grid version with cache
    makecache!(ch)
    for i in 1:numberofsteps
        last=deepcopy(ch.now)
        recexg(ch.now,ch.data,cach=ch.cache)
        if i%thinning==0
            push!(ch.history,last)
        end
        if ((~ silent) & ((i%1000)==0))
            println(i, "iterations complete")
        end
    end
end

function mcmc(ch::Fixedchain, numberofsteps=15000; thinning=1, silent=false) #All of the magic happens here. Fixed grid version
    ch.data=memorize(ch.data,ch.now.cop)
    for i in 1:numberofsteps
        last=deepcopy(ch.now)
        recexg(ch.now,ch.data)
        if i%thinning==0
            push!(ch.history,last)
        end
        if ((~ silent) & ((i%1000)==0))
            println(i, "iterations complete")
        end
    end
end


##
## ######## Some analysis ###########
##

function traceplot(ch::CopChain; burnin=1) #Make a traceplot of the chain
    energyhistory=[]
    for iter in ch.history[burnin:end]
        push!(energyhistory,iter.energy)
    end
    plot(energyhistory)
end


#graphing function for the mean. Works, but is very slow. Adjust 
function graphmean(ch::CopChain, burnin=4000, autocorr=100; logg=false)
    ch.now.cop.dim==2 || error("Plotting in more than two dimensions is not supported")
    alph=autocorr/length(ch.history[burnin:end])
    grid=zeros(100,100)
    for (i,iter) in enumerate(ch.history[burnin:end])
        ((i%autocorr)==0) && (grid+=forgraph(iter.cop,alpha=alph))
    end
    if logg
        grid=log.(grid)
    end
    heatmap(0:0.01:1,0:0.01:1,grid,c=cgrad([:black,:white]),legend=false, xticks=0:0.1:1, yticks=0:0.1:1)
end
