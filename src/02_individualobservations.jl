include("_common.jl")
using EcologicalNetworks
using StatsBase: mean, std
using Distributions: LogNormal, Categorical

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
            ict = sum(observations .== i )
            jct = sum(observations .== j)

            if ict > 0 && jct > 0
                tp += 1
            elseif ict > 0 || jct > 0
                fn += 1
            end
        end
    end
    return tp, fn
end

function getabundances(A)
    S = richness(A)
    # By trophic levels: 
    # Z = 2
    # abundances = Z.^[trophdict["s$i"]-1 for i in 1:S]
    # abundance_dist = abundances ./ sum(abundances)

    # By lognormal dist:
    abundances = rand(LogNormal(),S)
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
        sds[s] = mean(thesefrs);
        means[s] = std(thesefrs);
    end

    means, sds
end



buildplot(flexiblelinksmodel(10))

