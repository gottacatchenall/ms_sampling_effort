using EcologicalNetworks: degree_out, degree_in
using EcologicalNetworks
using Distributions
using DataFrames
using StatsBase
using Plots
using SparseArrays

function observe(A::T; falsepositive=0.1, falsenegative=0.3) where {T <: BipartiteNetwork}
    mat = Matrix(A.edges)
    for (i,el) in enumerate(mat)
        if el == 1 && rand() < falsenegative
            mat[i] = 0 
        elseif el == 0 && rand() < falsepositive
            mat[i] = 1
        end
    end
    return T(sparse(mat), A.T, A.B)
end

function observe(A::T; falsepositive=0.1, falsenegative=0.3) where {T <: UnipartiteNetwork}
    mat = Matrix(A.edges)
    for (i,el) in enumerate(mat)
        if el == 1 && rand() < falsenegative
            mat[i] = 0 
        elseif el == 0 && rand() < falsepositive
            mat[i] = 1
        end
    end
    return UnipartiteNetwork(mat)
end

function generate(truenet; richness=30, connectance=0.3, forbidden=nothing, kw...)
    O = observe(truenet; kw...)

    return truenet,O
end

function meandegcentrallity(net) 
    mean(values(centrality_degree(net)))
end

function entropy(net)
    return EcologicalNetworks.entropy(net)
end

function samplenetworks(A, fpr, fnr, property; numreplicates = 200)
    realstat = zeros(numreplicates)
    obstat = zeros(numreplicates)
    for r in 1:numreplicates
        metaweb = typeof(A) <: Function ? A(1) : A             

        real, obs = generate(metaweb, falsepositive=fpr,falsenegative=fnr)    
        realstat[r] = property(real)
        obstat[r] = property(obs)
    end
    
    error = sqrt.((realstat .- obstat).^2)
    return error
end

function get_error(netlist, fp, property)
    fnr = []
    meanerr = []
    err1σ = []
    fn_continuum = 0.0:0.02:1
     
    @showprogress for fn in fn_continuum
        thiserrs = []
        thissigmas = []
        for net in netlist
            error = samplenetworks(net, fp, fn, property)
            push!(thiserrs, mean(error))
            push!(thissigmas, sqrt(var(error)))
        end
        push!(fnr, fn)
        push!(meanerr, mean(thiserrs))
        push!(err1σ, mean(thissigmas))
    end
    return fnr, meanerr, err1σ
end
#=

    Running on MANGAL NETWORKS SECITON 

=#
include("get_mangal_data.jl")

fw, para, mutu, misc = mangaldata()

fpr = 0
fnr, para_err, para_sigma = get_error(para, fpr, η)
fnr, mutu_err, mutu_sigma = get_error(mutu, fpr, η)


plot(fnr, para_err)


for (mn,sg) in zip(mutu_errs, mutu_sigmas)
    plot!(pl, fn_continuum, mn, c=:red)
end
for (mn,sg) in zip(para_errs, para_sigmas)
    plot!(pl, fn_continuum, mn, c=:blue)
end
pl


fnr, para_err, para_sigma = get_error(para, fpr, connectance)
fnr, mutu_err, mutu_sigma = get_error(mutu, fpr, connectance)



scatter!(fnr, para_err)