using EcologicalNetworks
using Measures
using DataFrames
using Plots

# overload some functions 
Base.in(spname::String,net::NT) where {NT <: UnipartiteNetwork } = spname in net.S
Base.in(spname::String, net::NT) where {NT <: BipartiteNetwork } = spname in net.T || spname in net.B

function checkproductunipartite(M, marginalcount, jointcount)
    metaweb = adjacency(M)
    diffs = []
    joints = []
    margs = []

    S = richness(M)

    jointsum = sum(jointcount)
    margsum = sum(marginalcount)

    for i in 1:S
        P_i = marginalcount[i]/margsum
        for j in 1:S
            if i > j && metaweb[i,j] == 1
                P_j = marginalcount[j]/ margsum
                P_ij = jointcount[i,j] / jointsum # joint probab of seeing i and j
                thisdiff = P_ij - P_i*P_j
                push!(diffs, thisdiff)
                push!(joints, P_ij)
                push!(margs, P_i*P_j)
            end 
        end
    end
    return diffs, joints, margs
end

function checkproductbipartite(M, hostcounts, paracounts, jointcount)
    metaweb = adjacency(M)
    diffs = []
    joints = []
    margs = []
    numpara, numhost = size(metaweb)

    for i in 1:numpara
        P_i = paracounts[i]/sum(paracounts)
        for j in 1:numhost
            if metaweb[i,j] == 1
                P_j = hostcounts[j] / sum(hostcounts)
                P_ij = jointcount[i,j] / sum(jointcount)
                thisdiff = P_ij - P_i*P_j
                push!(diffs, thisdiff)
                push!(joints, P_ij)
                push!(margs, P_i*P_j)
            end 
        end
    end


    return diffs, joints, margs
end


function unipartitesplit(As)
    M = reduce(∪, As)
    speciespool = M.S
    S = richness(M)
    marginalcount = zeros(S)
    jointcount = zeros(S,S)
    for A in As 
        for (i,speciesi) in enumerate(speciespool)
            for (j, speciesj) in enumerate(speciespool)
                if i > j
                    marginalcount[i] +=  speciesi in A  
                    marginalcount[j] +=  speciesj in A
                    jointcount[i,j] += (speciesi in A && speciesj in A)
                end 
            end
        end
    end

    return checkproductunipartite(M,marginalcount, jointcount)
end

function bipartitesplit(Bs)
    M = reduce(∪, Bs)
    
    hostpool = M.B
    parapool = M.T
   
    host_marginalcount = zeros(length(hostpool))
    para_marginalcount = zeros(length(parapool))
    jointcount = zeros(length(parapool),length(hostpool), )

    for B in Bs 
        for (i,para) in enumerate(parapool)
            for (j,host) in enumerate(hostpool)
                    para_marginalcount[i] +=  para in B
                    host_marginalcount[j] +=  host in B  
                    jointcount[i,j] += (host in B && para in B)
            end
        end
    end

    return checkproductbipartite(M, host_marginalcount, para_marginalcount, jointcount)
end


function makejointmargplots(As)
    if typeof(As[begin]) <: UnipartiteNetwork
        margA, margB, jointAB = unipartitesplit(As);
    elseif typeof(As[begin]) <: BipartiteNetwork
        margA, margB, jointAB = bipartitesplit(As);
    end
end

As  = nz_stream_foodweb()
nz_diffs, nz_joints, nz_margs = makejointmargplots(As)

ids = map(i -> i.ID, filter(i -> contains(i.Reference, "Hadfield"), web_of_life()))
Bs = convert.(BipartiteNetwork, web_of_life.(ids))
hadfield_diffs, hadfield_joints, hadfield_margs = makejointmargplots(Bs)

using GLM

had_y = Float64.(hadfield_joints)
had_x = Float64.(hadfield_margs)
hadregression = lm(@formula(y~1+x), DataFrame(y=had_y, x = had_x))
had_α, had_β = coef(hadregression)

nz_y, nz_x = Float64.(nz_joints), Float64.(nz_margs)
nzregression = lm(@formula(y~1+x), DataFrame(y=nz_y, x = nz_x))
nz_α, nz_β, nz_σ = coef(nzregression)..., stderror(nzregression)[2]



fnt = font(20 , "Roboto")
hadfieldhist = histogram(hadfield_diffs,alpha=0.5,c=:dodgerblue4, xrotation=90,legend=:none, fontfamily=fnt, xlim=(-.00025, 0.0007), xticks=(-.00025:0.00025:0.0007),frame=:box)
ylabel!("frequency", tickfontfamily=fnt)
xlabel!("P(AB) - P(A)P(B)", tickfontfamily=fnt)
vline!([0], lw=2, c=:grey8, s=:dash)


hadfieldscat = scatter(hadfield_margs, hadfield_joints,
    fontfmaily=fnt,
    ma=0.05,
    ms=5, 
    mc=:dodgerblue4,
    msc=:dodgerblue4,
    legend=:none, 
    xlims = (0,0.0015),
    ylims = (0.0,0.0015),
    aspectratio=1.0, 
    frame=:box)

had_line_x, had_line_y = [i for i in 0:0.01:0.02], [had_α + i*had_β for i in 0:0.01:0.02]

plot!(had_line_x, had_line_y, lc=:dodgerblue4, lw=2, la=0.2)

plot!([0,1], [0,1], ls=:dash, la=9, lw=4, lc=:gray8)
ylabel!("P(AB)", tickfontfamily=fnt)
xlabel!("P(A)P(B)", tickfontfamily=fnt)



nzhist = histogram(nz_diffs,
    alpha=0.5, 
    c=:mediumpurple4, 
    xrotation=90,
    legend=:none, 
    fontfamily=fnt, 
    frame=:box)
ylabel!("frequency", tickfontfamily=fnt)
xlabel!("P(AB) - P(A)P(B)", tickfontfamily=fnt)
vline!([0], lw=2, c=:grey8, s=:dash)

nzscat = scatter(nz_margs, nz_joints,
    fontfmaily=fnt, 
    size=(400,400),
    xlims = (0,0.0002),
    ylims = (0,0.0002),
    ma=0.02,
    ms=5, 
    msw=1,
    mc=:mediumpurple4,
    msc=:gray40,
    legend=:none,
    frame=:box)


nz_line_x, nz_line_y = [i for i in  0:0.00001:0.0001], [nz_α + i*nz_β for i in  0:0.00001:0.0001]
nz_line_σ = nz_line_y .- [nz_α + i*(nz_β-nz_σ) for i in 0:0.00001:0.0001]


plot!(nz_line_x, nz_line_y, ribbon=nz_line_σ, lw=3, lc=:mediumpurple4, la=0.4)

plot!([0,1], [0,1], ls=:dash, la=9, lw=3, lc=:gray8)
ylabel!("P(AB)", tickfontfamily=fnt)
xlabel!("P(A)P(B)", tickfontfamily=fnt)


plot(hadfieldscat, hadfieldhist, nzscat, nzhist, layout = grid(2,2), margin=10mm, size=(800,800), dpi=300)


savefig("positiveassociations.png")
