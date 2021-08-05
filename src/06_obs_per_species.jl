using EcologicalNetworks
using Distributions



function samplingeffort_and_fnr(; 
    A = nichemodel(30, 0.1),
    numreplicates = 50,
    samplingeffort = vcat(1,25, 50:50:1500),
    λ=nothing
)

    sznet = typeof(A) <: Function ? A(1) : A             
    S = richness(sznet)
    sampling = []
    fnr_mean = []
    fnr_sigma = []

    for samp in samplingeffort
        thesefnrs = zeros(numreplicates)
        
        for rep in 1:numreplicates
            tp,fn = 0,0

            metaweb = typeof(A) <: Function ? A(1) : A             

            abundances = rand(LogNormal(),S)
            abundance_dist = abundances ./ sum(abundances)

            species_observered = zeros(Int64, samp)
            for i in 1:samp
                species_observered[i] = rand(Categorical(abundance_dist))
            end


            # check for A_ij = 1 where only one of i or j is observed
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


samp, fnr100_mean, fnr100_sd = samplingeffort_and_fnr(
    A=(_->nichemodel(100, 0.1)),
    numreplicates = 50)


samp, fnr50_mean, fnr50_sd = samplingeffort_and_fnr(
    A=(_->nichemodel(50, 0.1)),
    numreplicates = 50)


samp, fnr25_mean, fnr25_sd = samplingeffort_and_fnr(
    A=(_->nichemodel(25, 0.1)),
    numreplicates = 50)

fnt = font(20, "Roboto")

using ColorSchemes
cs2 = ColorScheme(range(colorant"dodgerblue", colorant"cyan4", length=3))

genplt = plot(dpi=300, la=0.3, ylims=(0,1),frame=:box, fontfamily=fnt,legendtitlefontsize=8, legendtitle="species richness", palette=cs2, colorbar=:right, size=(500,500))
plot!(genplt, samp, fnr50_mean, label="",c=cs2[1], ribbon=fnr50_sd, fa=0.3)
plot!(genplt, samp, fnr50_mean, label="", c=cs2[1], ribbon=fnr50_sd, fa=0.1)
scatter!(genplt, samp, fnr50_mean, ms=4, msw=2,label="50", msc=:dodgerblue, mc=:white)
plot!(genplt, samp, fnr25_mean, label="",c=cs2[2] ,ribbon=fnr25_sd, dpi=300, fa=0.3)
plot!(genplt, samp, fnr25_mean, label="",c=cs2[2] ,ribbon=fnr25_sd, dpi=300, fa=0.1)
scatter!(genplt, samp, fnr25_mean, ms=4,msw=2,label="25",msc=cs2[2], mc=:white)
plot!(genplt, samp, fnr100_mean, label="",ribbon=fnr100_sd, dpi=300, fa=0.3, c=cs2[3])
plot!(genplt, samp, fnr100_mean, label="",ribbon=fnr100_sd, dpi=300, fa=0.1, c=cs2[3])
scatter!(genplt, samp, fnr100_mean, ms=4,msw=2,label="100",msc=cs2[3], mc=:white)
yaxis!(genplt, "false negative rate")
xaxis!(genplt, "number of observations of individuals per species", xticks=0:200:1500, xlim=(0,1500))


