using GeoStats, JLD2
using ProgressMeter
#using GLMakie
using CairoMakie

include("utils/ghcn_data.jl")
include("utils/ols_trend.jl")

TREND_RANGE = 1.5e5
KRIG_RANGE = 1.5e5
DOY = 187
GRID = "1km"

MIN_TREND_STATS = 10
NLAGS = 10
MAX_NEIGH = 25
NODATA = -9999.

quants = [.05, .25, .5, .75, .95]
ints = [.99:-.01:.01...]

plot_labels = collect('A':'L')

f = Figure(size=(1500,1600))
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
    gamma = GeoStatsFunctions.fit(PentasphericalVariogram, vario, h->1)
    
    out_ras = Raster("data/ca_dem_1km.tif")
    x_out = [out_ras.xcoord[p[1]] for p in out_ras.data_locs]
    y_out = [out_ras.ycoord[p[2]] for p in out_ras.data_locs]
    out_pts = PointSet([(out_ras.xcoord[p[1]], out_ras.ycoord[p[2]])
                         for p in out_ras.data_locs])
    
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
    ax = Axis(g[1,1], aspect=1, xlabel="Easting (km)", ylabel="Northing (km)",
              title=string(plot_labels[i*3-2]), titlealign=:left)
    arr3 = plot_arr(out_ras, detrend_quants[:, 3])
    h = heatmap!(ax, out_ras.xcoord ./ 1e3, out_ras.ycoord ./ 1e3, arr3,
                 colormap=:inferno, colorrange=(-7.5, 47.0), rasterize=2)
    Colorbar(g[1,2], h, label="Tmax (deg C)")
    
    println(extrema(detrend_quants[:, 4] - detrend_quants[:, 2]))
    g = f[i,2] = GridLayout()
    ax = Axis(g[1,1], aspect=1, xlabel="Easting (km)", ylabel="Northing (km)",
              title=string(plot_labels[i*3-1]), titlealign=:left)
    arr1 = plot_arr(out_ras, detrend_quants[:, 4] - detrend_quants[:, 2])
    h = heatmap!(ax, out_ras.xcoord ./ 1e3, out_ras.ycoord ./ 1e3, arr1,
                 colormap=:turbo, colorrange=(1.1, 4.6), rasterize=2)
    Colorbar(g[1,2], h, label="Magnitude of 50% prediction interval (deg C)")
    #scatter!(ax, df.x, df.y, color=:black)
    
    println(extrema(detrend_quants[:, 5] - detrend_quants[:, 1]))
    g = f[i,3] = GridLayout()
    ax = Axis(g[1,1], aspect=1, xlabel="Easting (km)", ylabel="Northing (km)",
              title=string(plot_labels[i*3]), titlealign=:left)
    arr2 = plot_arr(out_ras, detrend_quants[:, 5] - detrend_quants[:, 1])
    h = heatmap!(ax, out_ras.xcoord ./ 1e3, out_ras.ycoord ./ 1e3, arr2,
                 colormap=:turbo, colorrange=(3.3, 10.5), rasterize=2)
    Colorbar(g[1,2], h, label="Magnitude of 90% prediction interval (deg C)")
end

save("figs/fig06.pdf", f)
#f
