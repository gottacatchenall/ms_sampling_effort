using EcologicalNetworks


Base.in(spname::String,net::NT) where {NT <: UnipartiteNetwork } = spname in net.S
Base.in(spname::String, net::NT) where {NT <: BipartiteNetwork } = spname in net.T || spname in net.B


# great function name
function dothing_bipartite(weblist, hostpool, parapool)

    host_marginalcount = zeros(length(hostpool))
    para_marginalcount = zeros(length(parapool))
    jointcount = zeros(length(hostpool), length(parapool))

    for B_site in weblist
        for (i,host) in enumerate(hostpool)
            for (j,para) in enumerate(parapool)
                if host in B_site
                    host_marginalcount[i] += 1
                end
                if para in B_site
                    para_marginalcount[j] += 1
                end
                if host in B_site && para in B_site
                    jointcount[i,j] += 1
                end
            end
        end
    end 

    metaweb = Matrix(M.edges)

    numpara, numhost = size(metaweb)

    diffs = []
    joints = []
    margs = []
    for i in 1:numpara
        P_i = para_marginalcount[i]/sum(para_marginalcount)
        for j in 1:numhost
            if metaweb[i,j] == 1
                P_j = host_marginalcount[j]/sum(host_marginalcount)
                P_ij = jointcount[j,i] / sum(jointcount)# joint probab of seeing i and j
                thisdiff = P_ij - P_i*P_j
                push!(diffs, thisdiff)
                push!(joints, P_ij)
                push!(margs, P_i*P_j)
            end 
        end
    end

    return diffs, joints, margs
end



using HypothesisTests
using Statistics
using Plots


ids = map(i -> i.ID, filter(i -> contains(i.Reference, "Hadfield"), web_of_life()))
B = convert.(BipartiteNetwork, web_of_life.(ids))
M = reduce(∪, B)

hostpool = M.B
parapool = M.T 
diffs, joints, margs = dothing_bipartite(B, M.B, M.T)


fnt = font(20, "Roboto")

hadfieldhist = histogram(diffs,alpha=0.5,c=:mediumpurple4, xrotation=90,legend=:none, fontfamily=fnt, xlim=(-.00025, 0.0007), xticks=(-.00025:0.00025:0.0007),frame=:box)
ylabel!("frequency", tickfontfamily=fnt)
xlabel!("P(AB) - P(A)P(B)", tickfontfamily=fnt)
vline!([0], lw=2, c=:grey8, s=:dash)


hadfieldscat = scatter(margs, joints,
    fontfmaily=fnt,
    ma=0.1,
    ms=3, 
    msw=2, 
    mc=:white,
    msc=:mediumpurple4,
    legend=:none, 
    xlims = (0,0.002),
    ylims = (0,0.002),
    aspectratio=1.0, 
    frame=:box)
plot!([0,1], [0,1], ls=:dash, la=9, lw=2, lc=:gray8)
ylabel!("P(AB)", tickfontfamily=fnt)
xlabel!("P(A)P(B)", tickfontfamily=fnt)

# ugh  
pvalue(OneSampleZTest(mean(diffs), std(diffs), length(diffs)))


####




# great function name
function dothing_unipartite(weblist)
    M = reduce(∪, weblist)
    speciespool = M.S
    S = length(speciespool)
    marginalcount = zeros(S)
    jointcount = zeros(S,S)

    for B_site in weblist
        for (i,speciesi) in enumerate(speciespool)
            for (j, speciesj) in enumerate(speciespool)
                if i > j
                    if speciesi in B_site
                        marginalcount[i] += 1
                    end
                    if speciesj in B_site
                        marginalcount[j] += 1
                    end
                    if speciesi in B_site && speciesj in B_site
                        jointcount[i,j] += 1
                    end
                end 
            end
        end
    end 

    metaweb = Matrix(M.edges)

    totalintct = sum(metaweb) 

    diffs = zeros(totalintct)
    joints = zeros(totalintct)
    margs = zeros(totalintct)
    

    jointsum = sum(jointcount)
    margsum = sum(marginalcount)

    cursor = 1
    for i in 1:S
        P_i = marginalcount[i]/margsum
        for j in 1:S
            if i > j && metaweb[i,j] == 1
                P_j = marginalcount[j]/ margsum
                P_ij = jointcount[i,j] / jointsum # joint probab of seeing i and j
                thisdiff = P_ij - P_i*P_j
                diffs[cursor] = thisdiff
                joints[cursor] = P_ij
                margs[cursor] = P_i*P_j
                cursor += 1
            end 
        end
    end

    return diffs, joints, margs
end




A  = nz_stream_foodweb()
nz_diffs, nz_joints, nz_margs = dothing_unipartite(A)


nzhist = histogram(nz_diffs,alpha=0.5, c=:crimson, xrotation=90,legend=:none, fontfamily=fnt, frame=:box)
ylabel!("frequency", tickfontfamily=fnt)
xlabel!("P(AB) - P(A)P(B)", tickfontfamily=fnt)
vline!([0], lw=2, c=:grey8, s=:dash)

nzscat = scatter(nz_margs, nz_joints,
    fontfmaily=fnt, 
    size=(400,400),
    xlims = (0,0.0002),
    ylims = (0,0.0002),
    ma=0.1,
    ms=3.5, 
    msw=1.25, 
    mc=:white,
    msc=:crimson,
    legend=:none,
    frame=:box)
plot!([0,1], [0,1], ls=:dash, la=9, lw=2, lc=:gray8)
ylabel!("P(AB)", tickfontfamily=fnt)
xlabel!("P(A)P(B)", tickfontfamily=fnt)



using Measures

plot(hadfieldscat,hadfieldhist, nzscat, nzhist, layout = grid(2,2), margin=10mm, size=(800,800), dpi=300)

savefig("positiveassociations.png")
