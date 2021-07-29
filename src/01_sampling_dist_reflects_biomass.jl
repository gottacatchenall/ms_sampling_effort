using DataFrames
using Distributions
using Plots
using EcologicalNetworks
using StatsBase


function samplingeffort_and_fnr(;
    A = nichemodel(S, 0.1),
    numreplicates = 50,
    samplingeffort = vcat(1,25, 50:50:1500),
    λ=nothing
)

    S = richness(A)
    sampling = []
    fnr_mean = []
    fnr_05 = []
    fnr_25 = []
    fnr_75 = []
    fnr_95 = []


    for samp in samplingeffort


        thesefnrs = zeros(numreplicates)
        for rep in 1:numreplicates
            tp,fn = 0,0

            abundances = rand(LogNormal(),S)
            abundance_dist = abundances ./ sum(abundances)

            species_observered = zeros(Int64, samp)
            for i in 1:samp
                species_observered[i] = rand(Categorical(abundance_dist))
            end

            for i in 1:S, j in 1:S
                if A[i,j] == 1
                    ict = sum(species_observered .== i )
                    jct = sum(species_observered .== j)

                    if ict > 0 && jct > 0
                        tp += 1
                    else 
                        fn += 1
                    end
                end
            end
            thesefnrs[rep] = fn/(fn+tp)
        end     

        push!(sampling, samp)

        push!(fnr_75, percentile(thesefnrs, 0.75))
        push!(fnr_25, percentile(thesefnrs, 0.25))
        push!(fnr_95, percentile(thesefnrs, 0.95))
        push!(fnr_05, percentile(thesefnrs, 0.05))

        push!(fnr_mean, mean(thesefnrs))
        
    end
    return sampling, fnr_mean, fnr_05, fnr_25, fnr_75, fnr_95
end

samp, fnr_mean, fnr_05, fnr_25, fnr_75, fnr_95 = samplingeffort_and_fnr(A=nichemodel(100, 0.1) ,numreplicates = 500)
plot(samp, fnr_mean, ribbon=(fnr_mean .- fnr_05, fnr_mean .-  fnr_95), dpi=300, fa=0.3, c=:dodgerblue, size=(700,500))
scatter!(samp, fnr_mean, ylim=(0,1), frame=:box,c=:white, ms=5, msw=2.5,msc=:dodgerblue, legend=:none, label="0.1",legendtitle="connectance")
yaxis!("false negative rate")
xaxis!("number of individual observations", xticks=0:100:1500, xlim=(0,1500))
savefig("samplingeffort_fnr.png")







#=

    DIFFERENT TOPOLOGICAL MODELS SECTION 

=#

samp, fnr_mean, fnr_05, fnr_25, fnr_75, fnr_95 = samplingeffort_and_fnr(A=nichemodel(100, 0.1) ,numreplicates = 500)



#=

    Running on MANGAL NETWORKS SECITON 

=#
include("get_mangal_data.jl")

fw, para, mutu, misc = mangaldata()


samps_per_fw = []
means_per_fw = []

@showprogress for thisfw in fw
    net = convert(UnipartiteNetwork, thisfw) 
    samp, fnr_mean, fnr_05, fnr_25, fnr_75, fnr_95 = samplingeffort_and_fnr(A=net,numreplicates = 500)
    push!(samps_per_fw, samp)
    push!(means_per_fw, fnr_mean)
end





plot(samp, fnr_mean, ribbon=(fnr_mean .- fnr_05, fnr_mean .-  fnr_95), dpi=300, fa=0.3, c=:dodgerblue, size=(700,500))
scatter!(samp, fnr_mean, ylim=(0,1), frame=:box,c=:white,msc=:dodgerblue, legend=:none, label="0.1",legendtitle="connectance")
yaxis!("false negative rate")
xaxis!("number of individual observations", xticks=0:100:1500, xlim=(0,1500))


convert(UnipartiteNetwork, fw[1])