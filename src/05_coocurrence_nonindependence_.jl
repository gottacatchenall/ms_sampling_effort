using EcologicalNetworks

ids = map(i -> i.ID, filter(i -> contains(i.Reference, "Hadfield"), web_of_life()))
B = convert.(BipartiteNetwork, web_of_life.(ids))
M = reduce(∪, B)

hostpool = M.B
parapool = M.T 

Base.in(spname::String,net::NT) where {NT <: UnipartiteNetwork } = spname in net.S
Base.in(spname::String, net::NT) where {NT <: BipartiteNetwork } = spname in net.T || spname in net.B

host_marginalcount = zeros(length(hostpool))
para_marginalcount = zeros(length(parapool))
jointcount = zeros(length(hostpool), length(parapool))


for B_site in B 
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

using Plots

heatmap(jointcount)

fullct = length(parapool)*length(hostpool)

diff = zeros(fullct)
cursor = 1
SCALING = 10^5
for (i,host) in enumerate(hostpool)
    P_i = host_marginalcount[i] / sum(host_marginalcount) # marginal probability of seeing i
    for (j,para) in enumerate(parapool)
        P_j = para_marginalcount[j]/sum(para_marginalcount)
        P_ij = jointcount[i,j] / sum(jointcount)# joint probab of seeing i and j
    end
end

using SparseArrays


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

using HypothesisTests
using Statistics
using Base: std

fnt = font(20, "Roboto")

hist = histogram(diffs,alpha=0.5, c=:teal, xrotation=90,legend=:none, fontfamily=fnt, xlim=(-.00025, 0.0007), xticks=(-.00025:0.00025:0.0007),frame=:box)
ylabel!("frequency", tickfontfamily=fnt)
xlabel!("P(AB) - P(A)P(B)", tickfontfamily=fnt, axisfontfamily="Computer Modern")
vline!([0], lw=2, c=:grey8, s=:dash)


scat = scatter(margs, joints,fontfmaily=fnt,alpha=0.1, ms=4, msw=2, mc=:mediumpurple4,msc=:mediumpurple4,legend=:none, xlim=(0,0.0015),ylim=(0,0.0015), aspectratio=1.0, frame=:box)
plot!([0,1], [0,1], ls=:dash, la=9, lw=2, lc=:gray8)
ylabel!("P(AB)", tickfontfamily=fnt)
xlabel!("P(A)P(B)", tickfontfamily=fnt, axisfontfamily="Computer Modern")

# ugh  
pvalue(OneSampleZTest(mean(diffs), std(diffs), length(diffs)))


using Measures

l = @layout [
    a{0.65w} b{0.9h}
]
plot(scat,hist, layout = l, margin=5mm, size=(600,400), dpi=300)

savefig("positiveassociations.png")