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
x = Float32.(copy(features[:, kept])) # Latent traits
y = Matrix(hcat(labels[kept])')       # Interaction bit


# ----------------------------------------------
# model setup

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
opt = ADAM() 


# This will setup the model for training,
# with 20% of data leftover for validation

function train(x, y, fnp; proportion=0.8,n_batches = 50000, batch_size=64, mat_at = 500    )
    epc = mat_at:mat_at:n_batches
    epc = vcat(1, epc...)

    trainlossvalue = zeros(Float64, n_batches)
    testlossvalue = zeros(Float64, n_batches)

    training_size = convert(Int64, floor(size(x, 2) * proportion))
   
    train = sort(sample(1:size(x, 2), training_size; replace=false))
    test = filter(i -> !(i in train), 1:size(x, 2))

# todo only add to training set
    training_with_falsenegatives = add_falsenegatives(y[:, train], fnp)


    data = (x[:, train], training_with_falsenegatives)


    data_test = (x[:, test], y[:, test])

    # This is the main training loop
    @showprogress for i in 1:n_batches
        # We pick a random batch out of the training set
        ord = sample(train, batch_size; replace=false)
        data_batch = (x[:, ord], y[:, ord])
        # If the training batch is too unbalanced, we draw another one
        while sum(data_batch[2]) < ceil(Int64, 0.25batch_size)
            ord = sample(train, batch_size; replace=false)
            data_batch = (x[:, ord], y[:, ord])
        end
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
    # Total number of cases
    n = tp .+ fp .+ tn .+ fn

    tpr = tp ./ (tp .+ fn)
    fpr = fp ./ (fp .+ tn)
    tnr = tn ./ (tn .+ fp)
    fnr = fn ./ (fn .+ tp)
    acc = (tp .+ tn) ./ (n)
    racc = ((tn .+ fp) .* (tn .+ fn) .+ (fn .+ tp) .* (fp .+ tp)) ./ (n .* n)
    bacc = ((tp ./ (tp .+ fn)) .+ (tn ./ (fp .+ tn))) ./ 2.0
    J = (tp ./ (tp .+ fn)) + (tn ./ (tn .+ fp)) .- 1.0
    ppv = tp ./ (tp .+ fp)

    thresholds = range(0.0, 1.0; length=500)
    thr_index = last(findmax(J))
    thr_final = thresholds[thr_index]

    pr = vec(m(features))    # Raw predictions
    imp = pr .>= thr_final   # Larger than threshold
    nwi = imp .& .!labels    # New

    P = copy(M)
    new_interactions = DataFrame(; parasite=String[], host=String[], cooc=Bool[], p=Float64[])
    cursor = 0
    for i in eachindex(species(P; dims=1)), j in eachindex(species(P; dims=2))
        cursor += 1
        if imp[cursor]
            P[i, j] = true
            if !M[i, j]
                push!(
                    new_interactions,
                    (species(P; dims=1)[i], species(P; dims=2)[j], cooc[cursor], pr[cursor]),
                )
            end
        end
    end
    return (roc=(fpr, tpr), pr=(tpr,ppv))
end


# ----------------------------------------------------------------------------
# here we go


function add_falsenegatives(input, fnr)
    noisey = similar(input)
    for (i,v) in enumerate(input)
        if v == 1 
            noisey[i] = rand() < fnr ? 0 : 1
        else 
            noisey[i] = input[i]
        end
    end
    return noisey
end




rocReal, prReal = train(x,basetruth, 0)
roc10, pr10 = train(x,labels_fnr10, 0.1)
roc25, pr25 = train(x,labels_fnr25, 0.25)
roc50, pr50 = train(x,labels_fnr50, 0.5)


cols = [:cyan2,:dodgerblue,:teal, :mediumpurple]

rocplt = plot(aspectratio=1, legend=:outerright, frame=:box, legendtitle="FNR", size=(750,500))
plot!(rocplt,rocReal[1], rocReal[2], lw=2.5, c=cols[1], label="")
scatter!(rocplt,rocReal[1], rocReal[2], lw=2.5, label="0", msw=1,msc=cols[1], mc=:white, ms=2)

plot!(rocplt,roc10[1], roc10[2], lw=2.5, c=cols[2], label="")
scatter!(rocplt,roc10[1], roc10[2], lw=2.5, c=cols[2], label="0.1", msw=1,msc=cols[2], mc=:white, ms=2)

plot!(rocplt,roc25[1], roc25[2], lw=2.5, c=cols[3], label="")
scatter!(rocplt,roc25[1], roc25[2], lw=2.5, label="0.25", msw=1,msc=cols[3], mc=:white, ms=2)

plot!(rocplt,roc50[1], roc50[2], lw=2.5, c=cols[4], label="")
scatter!(rocplt,roc50[1], roc50[2], lw=2.5, label="0.5", msw=1,msc=cols[4], mc=:white, ms=2)


plot!(rocplt, [0,1], [0,1], c=:grey, ls=:dash, la=0.8, aspectratio=1, label="random")
xaxis!(rocplt,"False positive rate", (0, 1))
yaxis!(rocplt,"True positive rate", (0, 1))

savefig(rocplt, "roc.png")

prplt =  plot(aspectratio=1, legend=:none, frame=:box, size=(500,500))
plot!(prplt,prReal[1], prReal[2], lw=3,c=cols[1])
scatter!(prplt,prReal[1], prReal[2], lw=2.5, label="0", msw=1,msc=cols[1], mc=:white, ms=2)

plot!(prplt,pr10[1], pr10[2], lw=3,c=cols[2])
scatter!(prplt,pr10[1], pr10[2], lw=2.5, label="0.1", msw=1,msc=cols[2], mc=:white, ms=2)

plot!(prplt,pr25[1], pr25[2], lw=3,c=cols[3])
scatter!(prplt,pr25[1], pr25[2], lw=2.5, label="0.25", msw=1,msc=cols[3], mc=:white, ms=2)

plot!(prplt,pr50[1], pr50[2], lw=3, c=cols[4])
scatter!(prplt,pr50[1], pr50[2], lw=2.5, label="0.5", msw=1,msc=cols[4], mc=:white, ms=2)


plot!(prplt, [0,1], [1,0], c=:grey, ls=:dash, la=0.5, aspectratio=1)
xaxis!(prplt,"True positive rate", (0, 1))
yaxis!(prplt,"Postive predictive value", (0, 1))

prplt


comb = plot(rocplt, prplt, size=(1250, 500), dpi=300)

savefig(comb, "rocpr_falsenegatives.png")

