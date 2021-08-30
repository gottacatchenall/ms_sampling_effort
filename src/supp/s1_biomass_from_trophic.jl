include("_common.jl")
using EcologicalNetworks
include("get_mangal_data.jl")

using StatsBase: mean, std
using Distributions: LogNormal, Categorical
using ColorSchemes, Colors
using Plots
using Measures

function observe(numobservations, relativeabundances)
    counts = zeros(length(relativeabundances))
    for i in 1:numobservations
        ind = rand(Categorical(relativeabundances))
        counts[ind] += 1
    end
    return counts
end

function countfp(A, observations)
    metaweb = adjacency(A)
    S = richness(A)
    tp, fn = 0,0
    for i in 1:S, j in 1:S
        if metaweb[i,j] == 1
            if observations[i] > 0 && observations[j] > 0
                tp += 1
            elseif observations[i] > 0 || observations[j] > 0
                fn += 1
            end
        end
    end
    return tp, fn
end

function getabundances(A)
    S = richness(A)
    Z = 2
    abundances = Z.^[trophdict["s$i"]-1 for i in 1:S]
    return  abundances ./ sum(abundances)
end

function runreplicate(numobservations, A)
    relativeabundances = getabundances(A)
    observations = observe(numobservations, relativeabundances)
    tp, fn = countfp(A, observations)

    fnr = fn/(fn+tp)
    return fnr
end

function buildplot(A;  
    numreplicates = 50,
    samplingeffort = vcat(1,25, 50:50:1500))


    means = zeros(length(samplingeffort))
    sds = zeros(length(samplingeffort))

    for (s, numobs) in enumerate(samplingeffort)
        thesefrs = zeros(numreplicates)
        for r in 1:numreplicates
            # sample metaweb if as is a function for generating,
            # otherwise it is an empirical web so just use it 
            metaweb = typeof(A) <: Function ? A(1) : A             
            thesefrs[r] = runreplicate(numobs, metaweb)
        end
        means[s] = mean(thesefrs);
        sds[s] = std(thesefrs);
    end

    means, sds
end


function makegeneratedplt(numreplicates)

    samp = vcat(1,25, 50:50:1500)
    mean250, sd250 = buildplot(flexiblelinksmodel(250), numreplicates=numreplicates)
    mean100, sd100 = buildplot(flexiblelinksmodel(100), numreplicates=numreplicates)
    mean50, sd50 = buildplot(flexiblelinksmodel(50), numreplicates=numreplicates)

    generatedplt = plot(
        dpi=300, 
        frame=:box, 
        fontfamily=fnt,
        ylims=(0,1), xlims=(0,1500))
    plot!(generatedplt, samp, mean50, label="", c=cs[1], ribbon=sd50, fa=0.3)
    plot!(generatedplt, samp, mean50, label="",  c=cs[1], ribbon=2sd50, fa=0.1)

    plot!(generatedplt, samp, mean100, label="", c=cs[2], ribbon=sd100, fa=0.3)
    plot!(generatedplt, samp, mean100, label="",  c=cs[2], ribbon=2sd100, fa=0.1)

    plot!(generatedplt, samp, mean250, label="", c=cs[3], ribbon=sd250, fa=0.3)
    plot!(generatedplt, samp, mean250, label="",  c=cs[3], ribbon=2sd250, fa=0.1)

    scatter!(generatedplt, samp, mean50, ms=4, msw=1, label="", msc=cs[1], mc=:white)
    scatter!(generatedplt, samp, mean100, ms=4, msw=1, label="", msc=cs[2], mc=:white)
    scatter!(generatedplt, samp, mean250, ms=4, msw=1, label="", msc=cs[3], mc=:white)
    title!(generatedplt,"a", titleloc=:left)
    yaxis!(generatedplt, "false negative rate", ylim=(0,1))
    xaxis!(generatedplt, "total number of individual observations", xticks=0:200:1500, xlim=(0,1500))

    return generatedplt
end

function makemangalplt(numreplicates)
    fw, para, mutu, misc = mangaldata()
    means_per_fw = []
    rich = [richness(f) for f in fw]

    @showprogress for thisfw in fw
        means, sds = buildplot(thisfw)
        push!(means_per_fw, means)
    end
    samp = vcat(1,25, 50:50:1500)

    means = means_per_fw[sortperm(rich)]

    cs1 = ColorScheme(range(colorant"dodgerblue", colorant"cyan4", length=length(rich)))
    empplt=plot(
        dpi=300, fontfamily=fnt, frame=:box, palette=cs1, legend=:none, colorbar=:right, size=(500,500))
    yaxis!(empplt, "false negative rate", ylim=(0,1))
    xaxis!(empplt, "total number of individual observations", xticks=0:200:1500, xlim=(0,1500))
    title!(empplt,"b", titleloc=:left)

    for mn in means
        plot!(empplt, samp, mn, la =0.3)
    end
    empplt
end


fnt = font(20, "Roboto")
cs = ColorScheme(range(colorant"dodgerblue", colorant"cyan4", length=3))

genplt = makegeneratedplt(200)
scatter!(genplt, [500], [500], ms=4, msw=1, label="50", msc=cs[1], mc=:white) # hack to get legend to work
scatter!(genplt, [500], [500], ms=4, msw=1, label="100", msc=cs[2], mc=:white) # hack to get legend to work
scatter!(genplt, [500], [500], ms=4, msw=1, label="250", msc=cs[3], mc=:white) # hack to get legend to work

empplt = makemangalplt(200)



emp = plot(yaxis=:none, ticks=:none, grid=:none, frame=:none)
plot(genplt, empplt, margin=5mm, dpi=300, size=(900,900))


savefig("figures/s1.png")