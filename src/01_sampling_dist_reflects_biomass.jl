using Distributions
using Plots
using EcologicalNetworks


S = 100 
A = nichemodel(S, 0.3)

numreplicates = 50
samplingeffort = [2^i for i in 2:10]

sampling = []
fnr = []

for samp in samplingeffort
    for rep in 1:numreplicates
        tp,fn = 0,0

        abundances = (rand(LogNormal(),S))
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
        # for each A_ij in the metaweb, check if i and j were observed
        tn = S^2 - links(A)

        push!(sampling, samp)
        push!(fnr,fn/(tn))
    end 
end

scatter(sampling, fnr, ylim=(0,1), frame=:box)
yaxis!("false negative rate")
xaxis!("number of individual observations")














histogram(species_observered, xlim=(1,100))
xlabel!("species by increasing abundance ->")
ylabel!("number of observations")