using Plots
using Distributions
using ColorSchemes
using Colors
numtimesobserved = 0:15

truenegativeprob_fnr10  = []
truenegativeprob_fnr25  = []
truenegativeprob_fnr50  = []
truenegativeprob_fnr75  = []
truenegativeprob_fnr90  = []


for t in numtimesobserved
    # its literally the birthday problem
    push!(truenegativeprob_fnr10, 1-(.1^t))
    push!(truenegativeprob_fnr25, 1-(0.25^t))
    push!(truenegativeprob_fnr50, 1-(0.5^t))
    push!(truenegativeprob_fnr75, 1-(0.75^t))
    push!(truenegativeprob_fnr90, 1-(0.9^t))
end


c1 = get(cs, 0.1)
c25 = get(cs, 0.25)
c5 = get(cs, 0.5)
c75 = get(cs, 0.75)
c9 = get(cs, 0.9)


function makeplot(cs)
    col(f) = get(cs,f)

    function add!(plt, x,y,sc)
        plot!(plt, x, y, lw=3, lc=col(sc),label="$(sc)")
        scatter!(plt, x, y, ms=5, msw=1.5, msc=col(sc), mc=:white, label="")    
    end
    plt = plot(legend=:outerright, legendtitle="false negative rate", xlim=(0, length(numtimesobserved)), ylim=(0,1), frame=:box, size=(800,500), dpi=300)
    xaxis!("Number of observed negatives")
    yaxis!("Probability interaction is a true negative")
    
    xticks!(0:15)
    add!(plt, numtimesobserved, truenegativeprob_fnr10, 0.1)
    add!(plt, numtimesobserved, truenegativeprob_fnr25, 0.25)
    add!(plt, numtimesobserved, truenegativeprob_fnr50, 0.5)
    add!(plt, numtimesobserved, truenegativeprob_fnr75, 0.75)
    add!(plt, numtimesobserved, truenegativeprob_fnr90, 0.9)
    return plt
end



pl = makeplot(colorschemes[:Paired_5])
pl

savefig("true_negatives.png")