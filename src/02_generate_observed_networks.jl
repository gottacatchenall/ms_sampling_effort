using EcologicalNetworks: degree_out, degree_in
using EcologicalNetworks
using Distributions
using DataFrames
using Plots


function observe(A; falsepositive=0.1, falsenegative=0.3)
    mat = Matrix(A.edges)
    for (i,el) in enumerate(mat)
        if el == 1 && rand() < falsenegative
            mat[i] = 0 
        elseif el == 0 && rand() < falsepositive
            mat[i] = 1
        end
    end
    return UnipartiteNetwork(mat)
end


function generate(; richness=30, connectance=0.3, forbidden=nothing, kw...)
    A = nichemodel(richness, connectance)
    O = observe(A; kw...)

    return A,O
end


real, obs = generate()

θ = (fp=[0, 0.01, 0.05,0.1], fn=0:0.01:0.9)

numreplicates = 500

totalct = length(θ[:fp]) * length(θ[:fn]) * numreplicates

falsepositiverate = zeros(totalct)
falsenegativerate = zeros(totalct)

real_connectance = zeros(totalct)
real_meanindegree = zeros(totalct)
real_meanoutdegree = zeros(totalct)

observed_connectance = zeros(totalct)
observed_meanoutdegree = zeros(totalct)
observed_meanindegree = zeros(totalct)


cursor = 1
for fp in θ[:fp]
    for fn in θ[:fn]
        for r in 1:numreplicates
            real, obs = generate(falsepositive=fp,falsenegative=fn)    
            falsepositiverate[cursor] = fp
            falsenegativerate[cursor] = fn

            real_connectance[cursor] = connectance(real)
            real_meanindegree[cursor] = mean(values(degree_out(real)))
            real_meanoutdegree[cursor] = mean(values(degree_in(real)))

            observed_connectance[cursor] = connectance(obs)
            observed_meanoutdegree[cursor] = mean(values(degree_out(obs)))
            observed_meanindegree[cursor] = mean(values(degree_in(obs)))
            cursor += 1
        end
    end
end

mse(x,y) = sqrt(sum((x .- y).^2))

mse_indegree = mse.(real_meanindegree,observed_meanindegree)
mse_outdegree = mse.(real_meanoutdegree, observed_meanoutdegree)
mse_connectance = mse.(real_connectance, observed_connectance)

df = DataFrame(fpr=falsepositiverate, fnr=falsenegativerate,mse_indeg=mse_indegree,mse_outdeg=mse_outdegree, mse_connectance=mse_connectance)



for fp in unique(df[!,:fpr])
    thisdf = filter(r -> r.fpr == fp && r.fnr == fn, df)
end

plot(plotlist...)


heatmap(Matrix(real.edges))
heatmap(Matrix(obs.edges))


fp0 = filter(r -> r.fpr == 0, df)

plot(frame=:box, xlim=(0,1), ylim=(0,1), aspectratio=1)
scatter!(fp0.fnr, fp0.mse_connectance, ma=0.01)