using JLD2, KernelDensity, Distributions
using CairoMakie

include("utils/result.jl")

quants = [0.005:0.005:0.995...]

res = load("output/res_187.jld2")["res"]
# pick some station with data
loc = res[200]
pred_quants = loc.quants

# two steps in quant vec == 0.01
cent_diff(ind, vec) = (0.01)/(2*(vec[ind]-vec[ind-1]))

function loc_slope(vec, win::Int=5)
    buff = win ÷ 2
    slopes = ones(length(vec) - 2*buff)
    for ind = buff+1:length(vec)-buff
        X = hcat(pred_quants[ind-buff:ind+buff], ones(win))
        y = quants[ind-buff:ind+buff]
        β̂ = X\y
        slopes[ind-buff] = β̂[1]
    end
    slopes
end

k = kde(pred_quants ; kernel=Epanechnikov)

high_lim = 0.28
f = Figure()
ax = Axis(f[1,1], xlabel="Tmax [C]", ylabel="Density",
          xgridvisible=false, ygridvisible=false)
ylims!(ax, 0, high_lim)
xlims!(ax, pred_quants[1], pred_quants[end])

cm = Reverse(:turbo)
mid = 100
for p = 5:5:99
    low = pred_quants[mid-p]
    high = pred_quants[mid+p]
    in_int = loc.true_val > low && loc.true_val < high
    α = in_int ? 1 : .25
    #α = .75
    ls = in_int ? :dash : :dot
    ls = in_int ? :solid : :solid
    col = in_int ? :blue : :red
    lines!(ax, [low, low], [0, high_lim], color=:blue, colormap=cm,
           colorrange=[1,99], alpha=α, linestyle=:solid, linewidth=1)
    lines!(ax, [high, high], [0, high_lim], color=:red, colormap=cm,
           colorrange=[1,99], alpha=α, linestyle=:solid, linewidth=1)
    println(p, " ", in_int)
end

lines!(ax, pred_quants, pdf(k, pred_quants), color=:black, linewidth=3, linestyle=:dot,
       label="Predictive PDF")
lines!(ax, [loc.true_val, loc.true_val], [0, high_lim], color=:black,
      linewidth=2.5, linestyle=:dash, label="Measured Tmax")
lines!(ax, [0,0], [0, 1], color=:blue, label="Interval lower bound")
lines!(ax, [0,0], [0, 1], color=:red, label="Interval upper bound")
axislegend(ax)
save("figs/fig03.pdf", f)
#f
