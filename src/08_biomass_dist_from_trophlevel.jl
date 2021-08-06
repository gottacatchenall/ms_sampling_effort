using DataFrames
using Distributions
using Plots
using EcologicalNetworks
using StatsBase


function samplingeffort_and_fnr(; 
    A = nichemodel(30, 0.1),
    numreplicates = 50,
    Z = 2,
    samplingeffort = vcat(1,25, 50:50:1500))

    sznet = typeof(A) <: Function ? A(1) : A             
    @show sznet
    S = richness(sznet)
    sampling = []
    fnr_mean = []
    fnr_sigma = []

    for samp in samplingeffort
        thesefnrs = zeros(numreplicates)
        for rep in 1:numreplicates
            tp,fn = 0,0

            metaweb = typeof(A) <: Function ? A(1) : A             


            #### TODO: alternatively sraw this from Z^(trophiclevl-1) as in yodzis
            ####

            # this shit is slow as hell for big networks 
            trophdict = trophic_level(metaweb)

            abundances = Z.^[trophdict["s$i"]-1 for i in 1:S]
            abundance_dist = abundances ./ sum(abundances)

            species_observered = zeros(Int64, samp)
            for i in 1:samp
                species_observered[i] = rand(Categorical(abundance_dist))
            end

            for i in 1:S, j in 1:S
                if metaweb[i,j] == 1
                    ict = sum(species_observered .== i )
                    jct = sum(species_observered .== j)

                    if ict > 0 && jct > 0
                        tp += 1
                    elseif ict > 0 || jct > 0
                        fn += 1
                    end
                end
            end
            thesefnrs[rep] = fn/(fn+tp)
        end     

        push!(sampling, samp)
        push!(fnr_sigma, sqrt(var(thesefnrs)))
        push!(fnr_mean, mean(thesefnrs))
        
    end
    return sampling, fnr_mean, fnr_sigma
end 

samp, fnr_mean, fnr_sd = samplingeffort_and_fnr(A=(_->nichemodel(50, 0.1)) ,numreplicates = 100)


fnt = font(20, "Roboto")

plot(samp, fnr_mean,fontfamily=fnt,  ribbon=fnr_sd, dpi=300, fa=0.3, c=:dodgerblue, size=(700,500))
scatter!(samp, fnr_mean, ylim=(0,1), frame=:box,c=:white, ms=5, msw=2.5,msc=:dodgerblue, legend=:none, label="0.1",legendtitle="connectance")
yaxis!("false negative rate")
xaxis!("number of individual observations", xticks=0:200:1500, xlim=(0,1500))
savefig("samplingeffort_fnr.png")




