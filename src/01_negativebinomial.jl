using Plots
using Distributions
using ColorSchemes
using Colors
using Measures

function add!(plt, fnr, col; maxtimesobserved=15)
    x,y = 1:maxtimesobserved, [1.0 - (fnr^t) for t in 1:maxtimesobserved]

    styl = (lw=3, ms=5, msw=1.5, mc=:white)

    plot!(plt, x, y, label="", lc=col; styl...)
    scatter!(plt, x, y, label=fnr, msc=col; styl...)    
end

numtimesobserved = 1:15

fnt = font(20, "Roboto")
cs = ColorSchemes.Paired_5

plt = plot(fonyfamily=fnt, 
    legendfontfamily=fnt,
    legend=:outerright, 
    xlim=numtimesobserved, 
    legendtitle="false negative rate", 
    ylim=(0,1), 
    frame=:box, size=(600,400), dpi=300, margin=10mm)
xaxis!("Number of observed negatives", fontfamily=fnt, xlims=(1,15), xticks=1:15)
yaxis!("Probability interaction is a true negative")

add!(plt, 0.1, cs[1])
add!(plt, 0.25, cs[2])
add!(plt, 0.5, cs[3])
add!(plt, 0.75, cs[4])
add!(plt, 0.9, cs[5])


savefig("figures/negativebinomial.png")
