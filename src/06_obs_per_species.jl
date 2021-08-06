using EcologicalNetworks
using Distributions
using ColorSchemes

function fig(numobs)
    relativeabund = []
    obs = []
    for relabd in vcat(0.001:0.0025:0.01, 0:0.01:0.2)
        push!(relativeabund, relabd)
        push!(obs, numobs/relabd)
    end
    relativeabund, obs
end

relabd, req1 = fig(1)
relabd, req10 = fig(10)
relabd, req100 = fig(100)
relabd, req500 = fig(500)
relabd, req1000 = fig(1000)

cs = ColorSchemes.tableau_sunset_sunrise

using Measures
ptsettings =  (mc=:white, msw=2.5, ms=5, lw=2.5)
fnt = font(20, "Roboto")
plot(frame=:box, fontfamily=fnt, legendtitlefontsize=9,size=(500,500), dpi=300,legend=:topright, margin=5mm, legendtitle="Goal # obs. of focal species")

plot!(relabd, req1, label="", lc=cs[1]; ptsettings...)
scatter!(relabd, req1, label="1", msc=cs[1]; ptsettings...)


plot!(relabd, req10, label="", lc=cs[2]; ptsettings...)
scatter!(relabd, req10, label="10", msc=cs[2]; ptsettings...)

plot!(relabd, req100, label="", lc=cs[3]; ptsettings...)
scatter!(relabd, req100, label="100", msc=cs[3]; ptsettings...)

plot!(relabd, req500, label="", lc=cs[4]; ptsettings...)
scatter!(relabd, req500, label="500", msc=cs[4]; ptsettings...)

plot!(relabd, req1000, label="", lc=cs[5]; ptsettings...)
scatter!(relabd, req1000, label="1000", msc=cs[5]; ptsettings...)

xaxis!("Relative abundance of focal species", xlim=(0,0.2), xticks=0:0.02:0.2)
yaxis!(:log10, "Expected needed observations of all species", yticks=[10^i for i in 0:5])


savefig("samplingeffort.png")