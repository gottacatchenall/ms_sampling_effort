using Distributions
using LinearAlgebra
using BioEnergeticFoodWebs

function draw_flexiblelinks_network(S)
    # MAP estimates from Macdonald et al 2020
    ϕ = 24.3
    μ = 0.086

    n = S^2 - (S-1)
    α = μ*ϕ
    β = (1-μ)*ϕ

    L = rand(BetaBinomial(n, α, β ))
    return BioEnergeticFoodWebs.nichemodel(S, L)
end




























p = model_parameters(net);
b = rand(size(net, 1));


# get dBdt 
BioEnergeticFoodWebs.dBdt(zeros(length(b)), b, p, 1)
 

# what i actually want here is the jacobian at the first step 
# applied to the same biomass vec with different adjacency matrices

s = simulate(p, b, start=0, stop=50)


function getjacobian(biomass, parameters)

    bm_matrix = zeros(eltype(biomass), (length(biomass), length(biomass)))
    rewire = (parameters[:rewire_method] == :ADBM) | (parameters[:rewire_method] == :Gilljam) | (parameters[:rewire_method] == :DS)
    costMat = rewire ? parameters[:costMat] : nothing
    BioEnergeticFoodWebs.fill_bm_matrix!(bm_matrix, biomass, parameters[:w], parameters[:A], parameters[:h], parameters[:functional_response], parameters[:ar]; rewire=rewire, costMat=costMat)
  
    # Available food
    F = zeros(eltype(biomass), (length(biomass), length(biomass)))
    BioEnergeticFoodWebs.fill_F_matrix!(F, bm_matrix, biomass, parameters[:Γh], parameters[:c], parameters[:ht], parameters[:functional_response])
    F
end


eigvals(getjacobian(b,p))



