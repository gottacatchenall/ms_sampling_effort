using DataFrames
using Distributions
using Plots
using EcologicalNetworks
using StatsBase


function samplingeffort_and_fnr(; 
    A = nichemodel(30, 0.1),
    numreplicates = 50,
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
            abundances = rand(LogNormal(),S)
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


### TODO: sclaing for connectance in these draws based on S-L relationship Macdonald et al 2020
samp, fnr_mean, fnr_sd = samplingeffort_and_fnr(A=(_->nichemodel(100, 0.1)) ,numreplicates = 300)

fnt = font(20, "Roboto")

plot(samp, fnr_mean,fontfamily=fnt,  ribbon=fnr_sd, dpi=300, fa=0.3, c=:dodgerblue, size=(700,500))
scatter!(samp, fnr_mean, ylim=(0,1), frame=:box,c=:white, ms=5, msw=2.5,msc=:dodgerblue, legend=:none, label="0.1",legendtitle="connectance")
yaxis!("false negative rate")
xaxis!("number of individual observations", xticks=0:200:1500, xlim=(0,1500))
savefig("samplingeffort_fnr.png")







#=

    Running on MANGAL NETWORKS SECITON 

=#
include("get_mangal_data.jl")

fw, para, mutu, misc = mangaldata()


samp, generated250_fnr_mean, generated250_fnr_sd = samplingeffort_and_fnr(A=(_->nichemodel(250, 0.1)),numreplicates = 500)
samp, generated100_fnr_mean, generated100_fnr_sd = samplingeffort_and_fnr(A=(_->nichemodel(100, 0.1)),numreplicates = 500)
samp, generated50_fnr_mean, generated50_fnr_sd = samplingeffort_and_fnr(A=(_->nichemodel(50, 0.1)), numreplicates = 500)
samp, generated25_fnr_mean, generated25_fnr_sd = samplingeffort_and_fnr(A=(_->nichemodel(25, 0.1)),numreplicates = 500)


plot(dpi=300, fontfamily=fnt, frame=:box, legendtitle="species")
plot!(samp, generated25_fnr_mean, ribbon=generated25_fnr_sd, label="25", fa=0.3)
plot!(samp, generated50_fnr_mean, ribbon=generated50_fnr_sd, label="50", fa=0.3)
plot!(samp, generated100_fnr_mean, ribbon=generated100_fnr_sd, label="100", fa=0.3)
plot!(samp, generated250_fnr_mean, ribbon=generated250_fnr_sd, label="250", fa=0.3)
yaxis!("false negative rate")
xaxis!("number of individual observations", xticks=0:200:1500, xlim=(0,1500))





samps_per_fw = []
means_per_fw = []

@showprogress for thisfw in fw
    samp, fnr_mean, fnr_sd = samplingeffort_and_fnr(A=thisfw,numreplicates = 200)
    push!(samps_per_fw, samp)
    push!(means_per_fw, fnr_mean)
end



using ColorSchemes

rich = [richness(f) for f in fw]
means = means_per_fw[sortperm(rich)]
cs1 = ColorScheme(range(colorant"dodgerblue", colorant"cyan4", length=length(rich)))
empplt=plot(samps_per_fw,means,dpi=300, fontfamily=fnt, la=0.3,frame=:box, palette=cs1, legend=:none, colorbar=:right, size=(500,500))
yaxis!(empplt, "false negative rate", ylim=(0,1))
xaxis!(empplt, "number of individual observations", xticks=0:200:1500, xlim=(0,1500))



fnt = font(20, "Roboto")


cs2 = ColorScheme(range(colorant"dodgerblue", colorant"cyan4", length=3))

genplt = plot(dpi=300, la=0.3, ylims=(0,1),frame=:box, fontfamily=fnt,legendtitlefontsize=8, legendtitle="species richness", palette=cs2, colorbar=:right, size=(500,500))
plot!(genplt, samp, generated50_fnr_mean, label="",c=cs2[1], ribbon=generated50_fnr_sd, fa=0.3)
plot!(genplt, samp, generated50_fnr_mean, label="", c=cs2[1], ribbon=2generated50_fnr_sd, fa=0.1)
scatter!(genplt, samp, generated50_fnr_mean, ms=4, msw=2,label="50", msc=:dodgerblue, mc=:white)
plot!(genplt, samp, generated100_fnr_mean, label="",c=cs2[2] ,ribbon=generated100_fnr_sd, dpi=300, fa=0.3)
plot!(genplt, samp, generated100_fnr_mean, label="",c=cs2[2] ,ribbon=2generated100_fnr_sd, dpi=300, fa=0.1)
scatter!(genplt, samp, generated100_fnr_mean, ms=4,msw=2,label="100",msc=cs2[2], mc=:white)
plot!(genplt, samp, generated250_fnr_mean, label="",ribbon=generated250_fnr_sd, dpi=300, fa=0.3, c=cs2[3])
plot!(genplt, samp, generated250_fnr_mean, label="",ribbon=2generated250_fnr_sd, dpi=300, fa=0.1, c=cs2[3])
scatter!(genplt, samp, generated250_fnr_mean, ms=4,msw=2,label="250",msc=cs2[3], mc=:white)
yaxis!(genplt, "false negative rate")
xaxis!(genplt, "number of individual observations", xticks=0:200:1500, xlim=(0,1500))


using Measures

l = @layout [grid(1,2)]
plot(genplt, empplt, dpi=300,size=(800,450), fontfamily=fnt,margin=4mm)
savefig("samplingdist.png")

yaxis!("false negative rate")
xaxis!("number of individual observations", xticks=0:200:1500, xlim=(0,1500))
savefig("samplingeffort_fnr_foodwebs.png")


