using Plots
using EcologicalNetworks
using Distributions

function draw_flexiblelinks_network(S)
    # MAP estimates from Macdonald et al 2020
    ϕ = 24.3
    μ = 0.086

    n = S^2 - (S-1)
    α = μ*ϕ
    β = (1-μ)*ϕ
    
    L = 0
    while L < 1 || (L > 0.5*S^2)
        L = rand(BetaBinomial(n, α, β))
    end
    return nichemodel(S, L)
end

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

function errormatrix(f;
    richnesses = 10:125,
    fnr = 0:0.01:1,
    nreps= 30
)

    mat = zeros(length(richnesses), length(fnr))

    for (si,S) in enumerate(richnesses)
        for (fni,p_fn) in enumerate(fnr)
            errorlist = zeros(nreps)

            for rep in 1:nreps

                truenetwork = draw_flexiblelinks_network(S)
                observednet = add_falsenegatives(truenetwork, p_fn)

                errorlist[rep] = sqrt((f(truenetwork) - f(observednet))^2)
            end
            mat[si, fni] = mean(errorlist)
        end
    end
    return mat 
end



connectance_err = errormatrix(connectance)
heatmap(connectance_err)
xaxis!("FNR")
yaxis!("richness")


centr_err = errormatrix(x-> mean(values(centrality_degree(x))))
heatmap(centr_err)

spectr_err = errormatrix(ρ)
heatmap(centr_err)