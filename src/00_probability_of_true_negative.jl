using Plots
using Distributions
using ColorSchemes
using Colors
using Measures
numtimesobserved = 0:15

fnt = font(20, "Computer Modern")

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



cs = ColorSchemes.Paired_5
plt = plot(fonyfamily=fnt, legendfontfamily=fnt,legend=:outerright, xlim=numtimesobserved, legendtitle="false negative rate", ylim=(0,1), frame=:box, size=(600,400), dpi=300)

function add!(plt, x,y,val, ct)
    plot!(plt, x, y, lw=3, lc=cs[ct],label="$(val)")
    scatter!(plt, x, y, ms=5, msw=1.5, msc=cs[ct], mc=:white, label="")    
end
xaxis!("Number of observed negatives", fontfamily=fnt)
yaxis!("Probability interaction is a true negative")

xticks!(0:15)
plot!(plt, xlims=(0,15))
add!(plt, numtimesobserved, truenegativeprob_fnr10, 0.1, 1)
add!(plt, numtimesobserved, truenegativeprob_fnr25, 0.25, 2)
add!(plt, numtimesobserved, truenegativeprob_fnr50, 0.5, 3)
add!(plt, numtimesobserved, truenegativeprob_fnr75, 0.75, 4)
add!(plt, numtimesobserved, truenegativeprob_fnr90, 0.9,5 )
plt


plot(plt, margin=5mm)
savefig("true_negatives.png")