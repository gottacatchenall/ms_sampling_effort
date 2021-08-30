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
degcentmat = errors(meandegcent, nreps=250)
specmat = errors(ρ, nreps=20) 

S = collect(30:5:150)
fnr = collect(0:0.02:1)

conplt = contourf(fnr, S, title="Connectance", colorbartitle="mse",connectmat, ylims=(30,50), xlims=(0,1), size=(400,400))
xlabel!(conplt,"false negative rate")
ylabel!(conplt,"species richness")



degplt = contourf(fnr, S, degcentmat, title="Degree Centrality", colorbartitle="mse", ylims=(30,50), xlims=(0,1), size=(400,400))
xlabel!(degplt,"false negative rate")
ylabel!(degplt,"species richness")

specplt = contourf(fnr, S, specmat, title="Spectral radius", colorbartitle="mse", ylims=(30,50), xlims=(0,1), size=(400,400))
xlabel!(specplt,"false negative rate")
ylabel!(specplt,"species richness")

emp = plot(axis=:none, ticks=:none, frame=:none)

using Measures

plot(conplt, degplt, specplt, emp, margin=5mm,layout=grid(2,2), dpi=300, size=(800,800))

savefig("s2.png")