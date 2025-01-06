using GeoStats, JLD2, KernelDensity
using Random
using CairoMakie
#using GLMakie

Random.seed!(1)

include("utils/ghcn_data.jl")
include("utils/ols_trend.jl")
include("utils/result.jl")

TREND_RANGE = 1.0e5
KRIG_RANGE = 1.0e5
DOY = 182
GRID = "1km"

MIN_TREND_STATS = 1
NLAGS = 12
MAX_NEIGH = 100
NODATA = -9999.

quants = [0.005:0.005:0.995...]

df, tmax = north_am(2022 ; grid=GRID)

# ids of two stations near Hancock, CA
inds = findall(in.(df.id, (["USC00043747", "USW00053119"],)))

# tmax for all stations on DOY
data = tmax[DOY, :]

# get rows and Tmax of stations
pred_stat = df[inds, :]
true_tmax = data[inds]

# drop stats with no data
has_data = (!ismissing).(data)
df_stats = df[has_data, :]
data = data[has_data]

# remove prediction locations
new_inds = findall(in.(df_stats.id, (pred_stat.id,)))
df_stats = df_stats[Not(new_inds), :]
data = data[Not(new_inds)]

# skip if prediction location in the same grid cell as another station with data
same_coord = innerjoin(pred_stat[:, [:gridx, :gridy]],
                       df_stats[:, [:gridx, :gridy]],
                       on=[:gridx, :gridy])
size(same_coord, 1) >= 1 && return NODATA

# detrend tmax, store in the df
df_stats.tmax = data
df_stats.dt_tmax, _ = ols_trend_domain(df_stats.x, df_stats.y,
                                       df_stats.elev, data ;
                                       MAXLAG=TREND_RANGE,
                                       MIN_TREND_STATS=MIN_TREND_STATS)

# normal score transform
stat_gr = georef(df_stats[:, [:x, :y, :dt_tmax]], (:x, :y))
norm_pipe = Quantile()
norm_stat_gr, norm_cache = apply(norm_pipe, stat_gr)

vario = EmpiricalVariogram(norm_stat_gr, :dt_tmax,
                           maxlag=KRIG_RANGE, nlags=NLAGS,)
gamma = GeoStatsFunctions.fit(PentasphericalVariogram, vario, h->1)

out_points = PointSet(GeoStats.Point.(pred_stat.gridx,
                                      pred_stat.gridy))

sol = norm_stat_gr |>
Interpolate(out_points, Kriging(gamma), prob=true,)

num_samp = 10000
out_points = PointSet(vcat(GeoStats.Point.(pred_stat.gridx,
                                           pred_stat.gridy),
                           GeoStats.Point.(df_stats.x,
                                           df_stats.y)))
lu = GeoStatsProcesses.LUMethod()
sim = rand(GeoStatsProcesses.GaussianProcess(gamma), out_points,
           norm_stat_gr, num_samp, lu)

trend = ols_trend_apply(pred_stat.gridx, pred_stat.gridy,
                        pred_stat.grid_elev,
                        df_stats.x, df_stats.y, df_stats.elev,
                        df_stats.tmax ; MAXLAG=TREND_RANGE,
                        MIN_TREND_STATS=MIN_TREND_STATS)
pred_stat.trend = trend

norm_quants = quantile.(sol.dt_tmax, (quants,))
norm_quants = reduce(hcat, norm_quants)'

function denorm_col(ind)
    gr = georef((dt_tmax=norm_quants[:, ind],), sol.geometry)
    revert(norm_pipe, gr, norm_cache).dt_tmax
end

denorm_quants = reduce(hcat, [denorm_col(i) for i in 1:length(quants)])
detrend_quants = denorm_quants .+ trend

krig_samp = hcat(rand.(sol.dt_tmax, num_samp)...)'

function denorm_samp(samp)
    gr = georef((dt_tmax=samp,), sol.geometry)
    revert(norm_pipe, gr, norm_cache).dt_tmax
end

denorm_krig_samp = reduce(hcat, [denorm_samp(krig_samp[:, i]) for i in 1:num_samp])
krig_samp = denorm_krig_samp .+ trend

# only first two points are simulated
denorm_sim_samp = reduce(hcat, [denorm_samp(sim[i].dt_tmax[1:2]) for i in 1:num_samp])
sim_samp = denorm_sim_samp .+ trend

bw = 1.
uk = kde(krig_samp' ; bandwidth=(bw,bw), npoints=(128, 128))
us = kde(sim_samp' ; bandwidth=(bw,bw), npoints=(128, 128))

low, high = extrema(vcat(krig_samp, sim_samp))
xs = LinRange(low, high, 100)
ys = LinRange(low, high, 100)

ik = InterpKDE(uk)
is = InterpKDE(us)

f = Figure(size=(600,590))
gtop = f[1, 1:2] = GridLayout()
rowsize!(gtop, 1, 130)
gbot = f[2, 1:2] = GridLayout()
rowsize!(gbot, 1, 320)
gright = f[1:2, 2] = GridLayout()
colsize!(gright, 1, 100)
gleft = f[1:2, 1] = GridLayout()
colsize!(gleft, 1, 320)

xlow = 30.75
xhigh = 44.25
ylow = 30.75
yhigh = 44.25

ax = Axis(f[1,1],  xticklabelsize=0, ylabel="Station 1 marginal density",
          ylabelsize=12, title="A", titlealign=:left)
lines!(ax, ik.itp.itp.ranges[1], sum(ik.itp.itp.itp, dims=1)[:] ./ sqrt(sum(ik.itp.itp.itp)) ; color=:blue)
lines!(ax, is.itp.itp.ranges[1], sum(is.itp.itp.itp, dims=1)[:] ./ sqrt(sum(is.itp.itp.itp)) ; color=:red)
xlims!(ax, xlow, xhigh)
ylims!(ax, 0., 0.22)

ax = Axis(f[2,2],  yticklabelsize=0, xlabel="Station 2 marginal density",
         xlabelsize=12, title="C", titlealign=:left)
lines!(ax, sum(ik.itp.itp.itp, dims=2)[:] ./ sqrt(sum(ik.itp.itp.itp)),
       ik.itp.itp.ranges[2] ; color=:blue, label="Kriging density")
lines!(ax, sum(is.itp.itp.itp, dims=2)[:] ./ sqrt(sum(is.itp.itp.itp)),
       is.itp.itp.ranges[2] ; color=:red, label="Simulation density")
ylims!(ax, ylow, yhigh)
xlims!(ax, 0., 0.22)

leg = Legend(f[1, 2], ax)
leg.framecolor = :white
ax = Axis(f[2,1], aspect=1, xlabel="Tmax station 1 [°C]",
          ylabel="Tmax station 2 [°C]", title="B", titlealign=:left)
contour!(ax, ik.itp.itp.ranges[1], ik.itp.itp.ranges[2], ik.itp.itp.itp ;
        levels=6, color=:blue)
contour!(ax, is.itp.itp.ranges[1], is.itp.itp.ranges[2], is.itp.itp.itp ;
         levels=6, color=:red)
xlims!(ax, xlow, xhigh)
ylims!(ax, ylow, yhigh)


k_cor = string(cor(krig_samp')[1, 2])[1:5]
s_cor = string(cor(sim_samp')[1, 2])[1:5]
text!(ax, 35.7, 31.5, ; text="Kriging corr: $k_cor", color=:blue, fontsize=12)
text!(ax, 35.7, 31, ; text="Simulation corr: $s_cor", color=:red, fontsize=12)

save("figs/fig07.pdf", f)

#f

mae(x, y) = mean(abs.(x .- y))

f = Figure(size=(400,300))
ax = Axis(f[1,1], ylabel="Density", xlabel="Tmax [°C]")
ylims!(ax, 0, .375)
stephist!(ax, sim_samp[2,:], normalization=:pdf, label="p(S1)")
println(mean(sim_samp[2,:]) - true_tmax[2], " ", std(sim_samp[2,:]))

good = abs.(sim_samp[1, :] .- true_tmax[1]) .< 2
stephist!(ax, sim_samp[2,good], normalization=:pdf, label="p(S1 | S2 error < 2)")
println(mean(sim_samp[2,good]) - true_tmax[2], " ", std(sim_samp[2,good]))

good = abs.(sim_samp[1, :] .- true_tmax[1]) .< 1
stephist!(ax, sim_samp[2,good], normalization=:pdf, label="p(S1 | S2 error < 1)")
println(mean(sim_samp[2,good]) - true_tmax[2], " ", std(sim_samp[2,good]))

lines!(ax, [true_tmax[2], true_tmax[2]], [0, .4], color=:black, linestyle=:dash, label="Measured S1 Tmax")
axislegend(ax, labelsize=8, patchsize=(12,8))

save("figs/fig08.pdf", f)
#f
