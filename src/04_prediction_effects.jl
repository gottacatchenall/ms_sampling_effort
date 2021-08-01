using Flux
using Random
using EcologicalNetworks
using StatsPlots
using MultivariateStats
using ProgressMeter
using BSON: @save, @load
using Statistics
using TSne
using DataFrames
using CSV: CSV
using ParallelKMeans
using JSON


# same cleaning as main.jl
for f in ["figures", "artifacts"]
    isdir(f) || mkdir(f)
end
Random.seed!(420);
theme(:bright)
ids = map(i -> i.ID, filter(i -> contains(i.Reference, "Hadfield"), web_of_life()))
B = convert.(BipartiteNetwork, web_of_life.(ids))
M = reduce(∪, B)
N = convert(UnipartiteNetwork, M)
K = EcologicalNetworks.mirror(N)
pc = fit(PPCA, Float64.(Array(K.edges)))
pr = MultivariateStats.transform(pc, Float64.(Array(K.edges)))

nf = 15
cooc = zeros(Bool, prod(size(M)))                  # H/P co-occurrence?
labels = zeros(Bool, prod(size(M)))                # H/P interaction?
features = zeros(Float64, (2 * nf, prod(size(M))))   # H/P latent traits
cursor = 0
for i in species(M; dims=1), j in species(M; dims=2)
    global cursor += 1
    # Interaction in the metaweb?
    labels[cursor] = M[i, j]
    # Values in the PCA space
    p_i = findfirst(i .== species(N))
    p_j = findfirst(j .== species(N))
    features[1:nf, cursor] .= pr[1:nf, p_i]
    features[(nf + 1):end, cursor] .= pr[1:nf, p_j]
    # Co-occurrence?
    for b in B
        if i in species(b)
            if j in species(b)
                cooc[cursor] = true
            end
        end
    end
end
kept = findall(cooc)
feat = Float32.(copy(features[:, kept])) # Latent traits
labs = Matrix(hcat(labels[kept])')       # Interaction bit


# ----------------------------------------------
# model setup


opt = ADAM() 


# This will setup the model for training,
# with 20% of data leftover for validation

function train(x, y, fnp; proportion=0.8,n_batches = 50000, batch_size=64, mat_at = 500    )
    m = Chain(
        # The first layer uses relu
        Dense(2nf, ceil(Int64, 2.4nf), relu),
        # The first dropout rate is 0.8 as we're using a lot of features
        Dropout(0.8),
        # All other layers are sigmoid with 0.6 dropout rate
        Dense(ceil(Int64, 2.4nf), ceil(Int64, 1.5nf), σ),
        Dropout(0.6),
        Dense(ceil(Int64, 1.5nf), ceil(Int64, 0.8nf), σ),
        Dropout(0.6),
        # The last layer has a single bit! P(parasite → host)
        Dense(ceil(Int64, 0.8nf), 1, σ),
    )
    loss(x, y) = Flux.mse(m(x), y)
    ps = Flux.params(m)

    epc = mat_at:mat_at:n_batches
    epc = vcat(1, epc...)

    trainlossvalue = zeros(Float64, n_batches)
    testlossvalue = zeros(Float64, n_batches)

    training_size = convert(Int64, floor(size(x, 2) * proportion))
   
    trainind = sort(sample(1:size(x, 2), training_size; replace=false))
    testind = filter(i -> !(i in trainind), 1:size(x, 2))

    newlabels = add_falsenegatives(y[:, trainind], fnp)
    data = (x[:, trainind], newlabels)

    data_test = (x[:, testind], y[:,testind])

    # This is the main training loop 
    for i in 1:n_batches
        # We pick a random batch out of the training set
        ord = sample(trainind, batch_size; replace=false)
        data_batch = (x[:, ord], y[:,ord])

        # This trains the model
        Flux.train!(loss, ps, [data_batch], opt)
        # We only save the loss at the correct interval
        if i in epc
            trainlossvalue[i] = loss(data...)
            testlossvalue[i] = loss(data_test...)
        end
    end

    # We get the predictions and observations for the testing dataset
    predictions = vec(m(data_test[1]))
    obs = vec(data_test[2])

    # And we pick thresholds in the [0,1] range
    thresholds = range(0.0, 1.0; length=500)

    # All this is going to be the components of the adjacency matrix at a given threshold
    tp = zeros(Float64, length(thresholds))
    fp = zeros(Float64, length(thresholds))
    tn = zeros(Float64, length(thresholds))
    fn = zeros(Float64, length(thresholds))

    # Main loop to get the four components
    for (i, thr) in enumerate(thresholds)
        pred = vec(predictions .>= thr)
        tp[i] = sum(pred .& obs)
        tn[i] = sum(.!(pred) .& (.!obs))
        fp[i] = sum(pred .& (.!obs))
        fn[i] = sum(.!(pred) .& obs)
    end

    tpr = tp ./ (tp .+ fn)
    fpr = fp ./ (fp .+ tn)
    ppv = tp ./ (tp .+ fp)


    return (roc=(fpr, tpr), pr=(tpr,ppv))
end


# ----------------------------------------------------------------------------
# here we go


function add_falsenegatives(input, fnr)
    noisey = zeros(size(input))
    for (i,v) in enumerate(input)
        if v == 1 
            noisey[i] = rand() < fnr ? 0 : 1
        else 
            noisey[i] = input[i]
        end
    end
    return vec(Bool.(noisey))
end




rocReal, prReal = [],[]
roc25, pr25 = [],[]
roc50, pr50 = [],[]
roc75, pr75 = [],[]

reps = 30
@showprogress for r in 1:reps 
    thisReal_roc, thisReal_pr = train(feat,labs, 0)
    push!(rocReal, thisReal_roc)
    push!(prReal, thisReal_pr)

    this25_roc, this25_pr = train(feat,labs, 0.25)
    push!(roc25, this25_roc)
    push!(pr25, this25_pr)
    
    this50_roc, this50_pr = train(feat,labs, 0.50)
    push!(roc50, this50_roc)
    push!(pr50, this50_pr)

    this75_roc, this75_pr = train(feat,labs, 0.75)
    push!(roc75, this75_roc)
    push!(pr75, this75_pr)
end 

using ColorSchemes
cols = [ColorSchemes.tableau_sunset_sunrise[i] for i in [1,2,3,4]]

fnt = font(20, "Roboto")


rocplt = plot()

for (thisx, thisy) in rocReal
    @show thisx, thisy
    plot!(rocplt, thisx, thisy, label="", lw=5,la=0.5, lc=cols[1])
end

for (thisx, thisy) in roc25
    plot!(rocplt, thisx, thisy, label="0.25", la=0.5, lc=cols[2])
end

for (thisx, thisy) in roc50
    plot!(rocplt, thisx, thisy, label="0.5",la=0.5, lc=cols[3])
end

for (thisx, thisy) in roc75
    plot!(rocplt, thisx, thisy, label="0.75", la=0.5, lc=cols[4])
end
rocplt















### OLD






rocplt = plot(aspectratio=1, fontfamily=fnt, legend=:outerright, frame=:box, legendtitle="FNR")
plot!(rocplt,rocReal[1], rocReal[2], lw=3, la=0.75, c=cols[1], label="0")
plot!(rocplt,roc10[1], roc10[2], lw=3, la=0.75,c=cols[2], label="0.1")
plot!(rocplt,roc25[1], roc25[2], lw=3, la=0.75,c=cols[3], label="0.25")
plot!(rocplt,roc50[1], roc50[2], lw=3, la=0.75, c=cols[4], label="0.5")


plot!(rocplt, [0,1], [0,1], c=:grey, ls=:dash, la=0.8, aspectratio=1, label="random")
xaxis!(rocplt,"False positive rate", (0, 1))
yaxis!(rocplt,"True positive rate", (0, 1))

savefig(rocplt, "roc.png")

prplt =  plot(aspectratio=1, fontfamily=fnt,legend=:outerright, frame=:box, legend_title="FNR")
plot!(prplt,prReal[1], prReal[2], lw=3,c=cols[1], la=0.8, label="0")
plot!(prplt,pr10[1], pr10[2], lw=3,c=cols[2], la=0.8, label="0.1")
plot!(prplt,pr25[1], pr25[2], lw=3,c=cols[3], la=0.8,label="0.25")
plot!(prplt,pr50[1], pr50[2], lw=3, c=cols[4], la=0.8, label="0.5")


plot!(prplt, [0,1], [1,0], c=:grey, ls=:dash, la=0.5, aspectratio=1, label="random")
xaxis!(prplt,"True positive rate", (0, 1))
yaxis!(prplt,"Postive predictive value", (0, 1))

prplt


comb = plot(rocplt, prplt, size=(900, 450), dpi=300, margin=3mm)

savefig(comb, "rocpr_falsenegatives.png")

