using GeoStats, JLD2
using ProgressMeter
#using GLMakie
using CairoMakie
using Proj
using Interpolations

include("utils/ghcn_data.jl")
include("utils/ols_trend.jl")

TREND_RANGE = 1.0e5
KRIG_RANGE = 1.0e5
GRID = "1km"

MIN_TREND_STATS = 10
NLAGS = 10
MAX_NEIGH = 100
NODATA = -9999.

quants = [.05, .25, .5, .75, .95]
ints = [.99:-.01:.01...]

plot_labels = collect('A':'L')

out_ras = Raster("data/ca_dem_1km.tif")
x_out = [out_ras.xcoord[p[1]] for p in out_ras.data_locs]
y_out = [out_ras.ycoord[p[2]] for p in out_ras.data_locs]
out_pts = PointSet([(out_ras.xcoord[p[1]], out_ras.ycoord[p[2]])
                     for p in out_ras.data_locs])
y, x = get_y_x_coords(out_ras.gt, out_ras.height, out_ras.width)
    
t = Proj.Transformation("EPSG:4326", "EPSG:3310", always_xy=true)
yrange = 32.5:.01:42.1
xrange = -124.8:.01:-113.5
g4326 = [(x, y) for y = yrange for x = xrange]
g3310 = [t(pt) for pt in g4326]
in_range(a) = a[1]>minimum(x) && a[1]<maximum(x) &&
              a[2]>minimum(y) && a[2]<maximum(y)
keep = in_range.(g3310)
g4326 = g4326[keep]
g3310 = g3310[keep]

f = Figure(size=(1500,1600))
f2 = Figure()
for (i, DOY) in enumerate([1, 91, 182, 274])
    df, tmax = north_am(2022 ; grid=GRID)
    df.tmax = tmax[DOY, :]
    df = df[(!ismissing).(df.tmax), :]
    
    df.dt_tmax, _, stat_counts = ols_trend_domain(df.x, df.y, df.elev, df.tmax
                                     ; MAXLAG=TREND_RANGE,
                                     MIN_TREND_STATS=MIN_TREND_STATS)
    
    # normal score transform
    stat_gr = georef(df[:, [:x, :y, :dt_tmax]], (:x, :y))
    norm_pipe = Quantile()
    norm_stat_gr, norm_cache = apply(norm_pipe, stat_gr)
    
    vario = EmpiricalVariogram(norm_stat_gr, :dt_tmax,
                               maxlag=KRIG_RANGE, nlags=NLAGS,)
    gamma = GeoStatsFunctions.fit(PentasphericalVariogram,
                                  vario, h->1)
    gamma_exp = GeoStatsFunctions.fit(ExponentialVariogram,
                                      vario, h->1)
    f2ax = Axis(f2[(i-1)÷2+1, (i+1)%2+1],
                xlabel="Distance [km]",
                ylabel="Semi-variogram",
                title=["A", "B", "C", "D"][i],
                titlealign=:left)
    scatter!(f2ax, getfield.(vario.abscissas, :val) ./ 1e3,
             vario.ordinates, color=:black, label="Empirical")
    lines!(f2ax, 1:100, gamma.(1:1000:100000), color=:blue,
          label="Pentaspherical")
    lines!(f2ax, 1:100, gamma_exp.(1:1000:100000), color=:red,
           label="Exponential")
    ylims!(f2ax, 0, 1.29)
    xlims!(f2ax, 0, 100)
    text!(f2ax, 15, 0.35 ; text="Nugget", color=:black)
    text!(f2ax, 15, 0.2 ;
          text=string(round(nugget(gamma); digits=2)),
          color=:blue)
    text!(f2ax, 15, 0.05 ;
          text=string(round(nugget(gamma_exp); digits=2)),
          color=:red)
    axislegend(f2ax, position=:rb, labelsize=8, markersize=1,
               patchsize=(10, 6))
    
    #pcoords = t.([(out_ras.xcoord[x[1]], out_ras.ycoord[x[2]])
    #              for x in out_ras.data_locs])
    #px = [x[1] for x in pcoords]
    #py = [x[2] for x in pcoords]
    sol = norm_stat_gr |>
          InterpolateNeighbors(out_pts, Kriging(gamma), prob=true,
                               neighborhood=MetricBall(KRIG_RANGE),
                               maxneighbors=MAX_NEIGH, minneighbors=1)
    
    if ismissing(sol.dt_tmax[1])
        println(ind, " ", doy)
    end
    
    trend = ols_trend_apply(x_out, y_out, out_ras.data,
                            df.x, df.y, df.elev,
                            df.tmax ; MAXLAG=TREND_RANGE,
                            MIN_TREND_STATS=MIN_TREND_STATS)
    
    norm_quants = quantile.(sol.dt_tmax, (quants,))
    norm_quants = reduce(hcat, norm_quants)'
    
    function denorm_col(ind)
        gr = georef((dt_tmax=norm_quants[:, ind],), sol.geometry)
        revert(norm_pipe, gr, norm_cache).dt_tmax
    end
    
    denorm_quants = reduce(hcat, [denorm_col(i) for i in 1:length(quants)])
    detrend_quants = denorm_quants .+ trend
    
    println(extrema(detrend_quants[:, 3]))
    g = f[i,1] = GridLayout()
    ax = Axis(g[1,1], aspect=1, xlabel="Longitude [°]",
              ylabel="Latitude [°]",
              title=string(plot_labels[i*3-2]), titlealign=:left)
    arr3 = plot_arr(out_ras, detrend_quants[:, 3])

    # Interpolate.jl requires increasing indices so some weird
    # indexing is required
    itp = interpolate((x, y[end:-1:1]), arr3[:, end:-1:1],
                      Gridded(Constant()))
    re = map(x->itp(x...), g3310)
    h = heatmap!(ax, [p[1] for p in g4326], [p[2] for p in g4326],
                 re, colormap=:inferno,
                 colorrange=(-10, 46.6), rasterize=true)
    Colorbar(g[1,2], h, label="Tmax [°C]")
    
    println(extrema(detrend_quants[:, 4] - detrend_quants[:, 2]))
    g = f[i,2] = GridLayout()
    ax = Axis(g[1,1], aspect=1, xlabel="Longitude [°]",
              ylabel="Latitude [°]",
              title=string(plot_labels[i*3-1]), titlealign=:left)
    arr1 = plot_arr(out_ras,
                    detrend_quants[:, 4]-detrend_quants[:, 2])
    itp = interpolate((x, y[end:-1:1]), arr1[:, end:-1:1],
                      Gridded(Constant()))
    re = map(x->itp(x...), g3310)
    h = heatmap!(ax, [p[1] for p in g4326], [p[2] for p in g4326],
                 re, colormap=:turbo, colorrange=(1.15, 4.55),
                 rasterize=true)
    Colorbar(g[1,2], h, label="Magnitude of 50% prediction interval [°C]")
    #scatter!(ax, df.x, df.y, color=:black)
    
    println(extrema(detrend_quants[:, 5] - detrend_quants[:, 1]))
    g = f[i,3] = GridLayout()
    ax = Axis(g[1,1], aspect=1, xlabel="Longitude [°]",
              ylabel="Latitude [°]",
              title=string(plot_labels[i*3]), titlealign=:left)
    arr2 = plot_arr(out_ras, detrend_quants[:, 5] - detrend_quants[:, 1])
    itp = interpolate((x, y[end:-1:1]), arr2[:, end:-1:1],
                      Gridded(Constant()))
    re = map(x->itp(x...), g3310)
    h = heatmap!(ax, [p[1] for p in g4326], [p[2] for p in g4326],
                 re, colormap=:turbo, colorrange=(3.6, 11.1),
                 rasterize=true)
    Colorbar(g[1,2], h, label="Magnitude of 90% prediction interval [°C]")
end

save("figs/fig06.pdf", f)
save("figs/appendix.pdf", f2)
#f
