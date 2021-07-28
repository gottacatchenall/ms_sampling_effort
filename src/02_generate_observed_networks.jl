using EcologicalNetworks: degree_out, degree_in
using EcologicalNetworks
using Distributions
using DataFrames
using StatsBase
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

function meanindeg(net) 
    mean(values(degree_in(net)))
end 

function meanoutdeg(net) 
    mean(values(degree_out(net)))
end

function meandegcentrallity(net) 
    mean(values(centrality_degree(net)))
end

function entropy(net)
    return EcologicalNetworks.entropy(net)
end

function samplenetworks(fpr, fnr, property; numreplicates = 500)
    realstat = zeros(numreplicates)
    obstat = zeros(numreplicates)
    for r in 1:numreplicates
        real, obs = generate(falsepositive=fpr,falsenegative=fnr)    
        realstat[r] = property(real)
        obstat[r] = property(obs)
    end
    
    error = sqrt.((realstat .- obstat).^2)
    return error
end


real, obs = generate()

function get_error(fp, property)

    fn_continuum = 0.0:0.02:1
    fnr = []
    meanerr = []
    err1σ = []
    err2σ = []
    for fn in fn_continuum
        error = samplenetworks(fp, fn, property)
        push!(fnr, fn)
        push!(meanerr, mean(error))

        σ = sqrt(var(error))

        push!(err1σ, σ)
        push!(err2σ, 2σ)
    end

    return fnr, meanerr, err1σ, err2σ
end


fnr, mean_connectance_err, connectance_1sigma, connectance_2sigma  = get_error(0.0, connectance)
connectplt = plot(frame=:box, legend=:none, size=(400,400))
plot!(connectplt, fnr, mean_connectance_err, ribbon=connectance_1sigma, fc=:dodgerblue,fa=0.5)
plot!(connectplt,fnr, mean_connectance_err, ribbon=connectance_2sigma, fc=:dodgerblue, fa=0.2)
scatter!(connectplt,fnr, meanerr_fp0, c=:white, ms=3, msw=1.5, msc=:dodgerblue, label="FPR = 0")
xaxis!(connectplt,"false negative rate", xlims=(0,1))
yaxis!(connectplt,"Error in estimated connectance", ylims=(0,1))


fnr, indeg_err, indeg_1sigma, indeg_2sigma  = get_error(0.0, meanindeg)
degplt = plot(frame=:box, legend=:none, size=(400,400))
plot!(degplt,fnr, indeg_err, ribbon=indeg_1sigma,lc=:seagreen4, fc=:seagreen4,fa=0.5)
plot!(degplt,fnr, indeg_err, ribbon=indeg_2sigma, lc=:seagreen4,fc=:seagreen4, fa=0.2)
scatter!(degplt,fnr, indeg_err, c=:white, ms=3, msw=1.5,msc=:seagreen4, label="FPR = 0")
xaxis!(degplt,"false negative rate", xlims=(0,1))
yaxis!(degplt,"Error in estimated degree", ylims=(0,10))



fnr, degcentrality_err,  degcentrality_1sigma,  degcentrality_2sigma  = get_error(0.0, meandegcentrallity)
degcentrplt = plot(frame=:box, legend=:none, size=(400,400))
plot!(degcentrplt, fnr, degcentrality_err, ribbon=degcentrality_1sigma,lc=:teal, fc=:teal,fa=0.5)
plot!(degcentrplt, fnr, degcentrality_err, ribbon=degcentrality_2sigma, lc=:teal,fc=:teal, fa=0.2)
scatter!(degcentrplt, fnr, degcentrality_err, c=:white, ms=3, msw=1.5,msc=:teal, label="FPR = 0")
xaxis!(degcentrplt,"false negative rate", xlims=(0,1))
yaxis!(degcentrplt,"Error in estimated mean degree centality", xlims=(0,1))




fnr, ent_err, ent_1sigma, ent_2sigma  = get_error(0.0, entropy)
entplt = plot(frame=:box, legend=:none, ylim=(0,7), size=(400,400))
plot!(entplt,fnr, ent_err, ribbon=ent_1sigma,lc=:mediumpurple4, fc=:mediumpurple4,fa=0.5)
plot!(entplt,fnr, ent_err, ribbon=ent_2sigma, lc=:mediumpurple4,fc=:mediumpurple4, fa=0.2)
scatter!(entplt,fnr, ent_err, c=:white, ms=3, msw=1.5,msc=:mediumpurple4, label="FPR = 0")
xaxis!(entplt,"false negative rate", xlims=(0,1))
yaxis!(entplt,"Error in estimated entropy")


plot(connectplt, degplt, degcentrplt, entplt, size=(700,700), dpi=300, padding=10)


savefig("properties_error.png")