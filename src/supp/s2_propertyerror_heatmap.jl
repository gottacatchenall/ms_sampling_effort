using Plots
using EcologicalNetworks
using StatsBase
using Statistics: std

include("_common.jl")

function add_falsenegatives(truenetwork, p_fn)
    mat = adjacency(truenetwork)
    for (i,el) in enumerate(mat)
        if el == 1 && rand() < p_fn
            mat[i] = 0 
        elseif el == 0 && rand() < p_fn
            mat[i] = 1
        end
    end
    return UnipartiteNetwork(mat)
end

function errors(f;
    S = 30:5:150,
    fnr = 0:0.02:1,
    nreps= 30
)

    err = zeros(length(S),length(fnr))

    for (si,s) in enumerate(S)
        for (fni,p_fn) in enumerate(fnr)
            errorlist = zeros(nreps)
            for rep in 1:nreps

                truenetwork = flexiblelinksmodel(s)
                observednet = add_falsenegatives(truenetwork, p_fn)
                errorlist[rep] = sqrt((f(truenetwork) - f(observednet) )^2)
            end
            err[si, fni] = mean(filter(x -> !isnan(x), errorlist))
        end

    end

    return err
end

meandegcent(x) = mean(collect(values(centrality_degree(x))))


connectmat= errors(connectance, nreps=100)
degcentmat = errors(meandegcent, nreps=100)
specmat = errors(ρ, nreps=20) 


heatmap(connectmat)

heatmap(degcentmat)