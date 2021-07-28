using Distributions
using Plots

S = 100 

abundances = sort(rand(LogNormal(),S))

abundance_dist = abundances ./ sum(abundances)

observations = 1000
species_observered =zeros(Int64, observations)
for i in 1:observations
    species_observered[i] = rand(Categorical(abundance_dist))
end




histogram(species_observered, xlim=(1,100))
xlabel!("increasing abundance ->")
ylabel!("number of observations")