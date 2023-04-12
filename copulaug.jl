##
## Grid continuous copulas, plus procedures to change them.
## Some procedures have extra parameters, aimed at making reversible-jump MCMC easier.
##

## If this says anything about baboons, you are trying to run something that isn't implemented yet

# using Distributions
using Plots

#
# Define the Copula structure
#

mutable struct Copula #Asymptotically approach any continuous copula
    dim::Int
    divs::Array{Array{Float64}}
    sizes::Array{Float64}
    probs::Array{Float64}
end
export Copula

#In order to do anything, the interval must besplit along 2 different dimensions, so this is part of the initialization procedure
function makefirsttwosplits(cop::Copula)
    #pick 2 dimensions
    dim1=rand(1:cop.dim)
    dim2=rand(1:(cop.dim-1))
    dim2>=dim1 && (dim2+=1)
    #pick the locations to split
    wheresplit1=rand()
    wheresplit2=rand()
    cop.divs[dim1]=hcat(wheresplit1, cop.divs[dim1])
    cop.divs[dim2]=hcat(wheresplit2, cop.divs[dim2])
    #Calculate the probabilities of the split cells
    dimofones=ones(Int,cop.dim) #To make sure we have the right number of dimensions
    cop.probs[dimofones...]=wheresplit1*wheresplit2
    cop.probs=cat(dims=dim1, cop.probs, [(1-wheresplit1)*wheresplit2])
    blump=[wheresplit1*(1-wheresplit2)]
    blump=cat(dims=dim1,blump,(1-wheresplit1)*(1-wheresplit2))
    cop.probs=cat(dims=dim2, cop.probs, blump)
    cop.sizes=copy(cop.probs)
end

function Copula(n::Int) #constructor function for copulas
    dim=n
    divs=[] #where we store the locations of all the divisions This is a list of lists, where the nth list is the divisions along the nth dimension
    for i in 1:n
        push!(divs,[1.0])
    end
    dimofones=ones(Int, dim)
    probs=ones(dimofones...)
    sizes=ones(dimofones...)
    cops=Copula(dim, divs, probs,sizes)
    makefirsttwosplits(cops) #initialize it with two splits already made
    return cops
end

#change the values in a rectangle
function rectangle(cop::Copula, dim1=0, dim2=0, loc1=0, newplace1=0, newplace2=0; rewrite=false)
    #pick 2 dimensions with at least one split each
    if dim1==0 #If dimensions were provided, ignore this part
        validdims=[]
        for i in 1:cop.dim
            (cop.divs[i]==[1.]) || (push!(validdims, i))
        end
        dim1=rand(1:length(validdims))
        dim2=rand(1:(length(validdims)-1))
        dim2>=dim1 && (dim2+=1)
        dim1=validdims[dim1]
        dim2=validdims[dim2]
    end
    #pick the cells to change from the array of probabilities
    if loc1==0 #If the locations were provided, ignore this
        sz=size(cop.probs)
        loc1=[]
        for s in sz
            push!(loc1,rand(1:s))
        end
    end
    indices1=CartesianIndex(loc1...)
    value1=cop.probs[indices1]
    loc2=Array(loc1)
    if newplace1==0 #If the locations were provided, ignore this
        newplace1=rand(1:length(cop.divs[dim1])-1)
        newplace1>=loc1[dim1] && (newplace1+=1)
    end
    loc2[dim1]=newplace1
    indices2=CartesianIndex(loc2...)
    value2=cop.probs[indices2]
    loc3=loc2
    if newplace2==0 #If the locations were provided, ignore this
        newplace2=rand(1:length(cop.divs[dim2])-1)
        newplace2>=loc1[dim2] && (newplace2+=1)
    end
    loc3[dim2]=newplace2
    indices3=CartesianIndex(loc3...)
    value3=cop.probs[indices3]
    loc4=loc1
    loc4[dim2]=newplace2
    indices4=CartesianIndex(loc4...)
    value4=cop.probs[indices4]
    #change them
    epsilon=rand()*(min(value1,value3)+min(value2,value4))+(-min(value2,value4)) #Uniform over the range
    value2+=epsilon
    value4+=epsilon
    value1-=epsilon
    value3-=epsilon
    probies=copy(cop.probs)
    probies[indices1]=value1
    probies[indices2]=value2
    probies[indices3]=value3
    probies[indices4]=value4
    rewrite && (cop.probs=probies)
    return(probies,dim1,dim2,indices1,indices2,indices3,indices4)
end

function forgraph(cop::Copula; alpha=1) #Works but is slow. Need to optimize
    cop.dim==2 || (error("BABOONS don't know how to plot in more than two dimensions!!!"))
    forgrid=ones(100,100)
    dens=cop.probs ./ cop.sizes
    for (i,hor) in enumerate(0.01:0.01:1)
        indexhor=findfirst(l->cop.divs[1][l]>=hor,1:100)
        for (j,vert) in enumerate(0.01:0.01:1)
            indexvert=findfirst(l->cop.divs[2][l]>=vert,1:100)
            forgrid[i,j]=dens[indexhor,indexvert]
        end
    end
    return forgrid/alpha
end

##Faster but throwing weird exceptions. Baboons for the time being.
#function forgraph(cop::Copula; alpha=1)
#    cop.dim==2 || (error("BABOONS don't know how to plot in more than two dimensions!!!"))
#    forgrid=[]
#    dens=cop.probs ./ cop.sizes
#    wherenowhor=1
#    graphstillessh=0
#    graphstillesshp=0
#    while cop.divs[1][wherenowhor]<1
#        blockhors=[]
#        while ((0:0.01:1)[graphstillessh])<wherenowhor
#            graphstillessh+=1
#        end
#        wherenowvert=1
#        graphstillessv=0
#        graphstillesvp=0
#        while cop.divs[2][wherenowvert]<=1
#            while ((0:0.01:1)[graphstillessv])<wherenowvert
#                graphstillessv+=1
#            end
#        end
#        blockverts=ones(graphstillessh-graphstillesshp,graphstillessv-graphstillesvp)*dens[wherenowhor,wherenowvert]
#        blockhors=vcat(blockhors,blockverts)
#    end
#    forgrid=hcat(forgrid,blockhors)
#    return forgrid*alpha
#end

function graph(cop::Copula)
    forg=forgraph(cop)
    heatmap(0:0.01:1,0:0.01:1,forg,c=cgrad([:black,:white]),legend=false, xticks=0:0.1:1, yticks=0:0.1:1)
end
