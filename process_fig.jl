using GeoIO, GeoStats
using CairoMakie

include("utils/ghcn_data.jl")
include("utils/ols_trend.jl")

MIN_TREND_STATS = 1
NLAGS = 12
DOY = 187

bound_fn = "data/ca_bounds/ca_bounds_3310.shp"
bound = GeoIO.load(bound_fn)

df, tmax = north_am(2022 ; grid="1km")
in_ca = check_in_bound(df.lon, df.lat, :ca)
df = df[in_ca, :]
tmax = tmax[DOY, in_ca]
df.tmax = tmax
df = df[(!ismissing).(df.tmax), :]
df.tmax = convert(Vector{Float64}, df.tmax)

wong = Makie.wong_colors()
bar_gap = 2
f = Figure(size=(1100, 340))
axis_label_size = 10
title_size = 10
tick_label_size = 9

g1 = f[1, 1] = GridLayout()
ax = Axis(g1[1,1], aspect=1, xticklabelsize=tick_label_size, yticklabelsize=tick_label_size,
          ylabel="Northing (km)", ylabelsize=axis_label_size,
          yticks=(-4e5:4e5:4e5, string.([-400:400:400...])),
          xticks=(-4e5:4e5:4e5, string.([-400:400:400...])),
          xlabel="Easting (km)", xlabelsize=axis_label_size,
          limits=(-4.2e5, 5.7e5, -6.6e5, 4.9e5),
          title="A", titlesize=title_size, titlealign=:left)

function plot_outline(ax)
    for geom in bound.geometry[1].geoms 
        lines!(ax, map(x -> (x.coords[1], x.coords[2]), geom.rings[1].vertices),
               color=:black, label="Study area", linewidth=.5)
    end
end

plot_outline(ax)
s = scatter!(ax, df.x, df.y, markersize=2, color=df.tmax, colormap=:inferno)
Colorbar(g1[1, 2], s, ticklabelsize=tick_label_size,)
colgap!(g1, bar_gap)

ax = Axis(f[2,1], yticklabelsize=tick_label_size, xticklabelsize=tick_label_size, ylabel="Density",
          ylabelsize=axis_label_size,
          xlabel="Tmax [C]", xlabelsize=axis_label_size,
          title="B", titlesize=title_size, titlealign=:left)
ylims!(ax, 0, 0.082)
hist!(f[2,1], df.tmax, normalization=:pdf, color=wong[1])


TREND_RANGE = 3.5e5
KRIG_RANGE = 3.5e5

df.dt_tmax, _ = ols_trend_domain(df.x, df.y, df.elev, df.tmax
                                 ; MAXLAG=TREND_RANGE,
                                   MIN_TREND_STATS=MIN_TREND_STATS)

g2 = f[1, 2] = GridLayout()
ax = Axis(g2[1,1], aspect=1, xticklabelsize=tick_label_size, yticklabelsize=tick_label_size,
          yticks=(-4e5:4e5:4e5, string.([-400:400:400...])),
          xticks=(-4e5:4e5:4e5, string.([-400:400:400...])),
          xlabel="Easting (km)", xlabelsize=axis_label_size,
          limits=(-4.2e5, 5.7e5, -6.6e5, 4.9e5),
          title="C", titlesize=title_size, titlealign=:left)
ax.yticklabelsize = 0
plot_outline(ax)
s = scatter!(ax, df.x, df.y, markersize=2, color=df.dt_tmax,
             colormap=:viridis)
Colorbar(g2[1, 2], s, ticklabelsize=tick_label_size,)
colgap!(g2, bar_gap)

ax = Axis(f[2,2], yticklabelsize=tick_label_size, xticklabelsize=tick_label_size, xlabel="De-trended Tmax [C]", xlabelsize=axis_label_size,
          title="D", titlesize=title_size, titlealign=:left)
ylims!(ax, 0, 0.22)
hist!(f[2,2], df.dt_tmax, normalization=:pdf, color=wong[2])


stat_gr = georef(df[:, [:x, :y, :dt_tmax]], (:x, :y))
norm_pipe = Quantile()
norm_stat_gr, norm_cache = apply(norm_pipe, stat_gr)

g3 = f[1, 3] = GridLayout()
ax = Axis(g3[1,1], aspect=1, xticklabelsize=tick_label_size, yticklabelsize=tick_label_size,
          yticks=(-4e5:4e5:4e5, string.([-400:400:400...])),
          xticks=(-4e5:4e5:4e5, string.([-400:400:400...])),
          xlabel="Easting (km)", xlabelsize=axis_label_size,
          limits=(-4.2e5, 5.7e5, -6.6e5, 4.9e5),
          title="E", titlesize=title_size, titlealign=:left)
ax.yticklabelsize = 0
plot_outline(ax)
pts_x = [coordinates(x)[1] for x in norm_stat_gr.geometry]
pts_y = [coordinates(x)[2] for x in norm_stat_gr.geometry]
s = scatter!(ax, pts_x, pts_y, markersize=2, color=norm_stat_gr.dt_tmax,
             colormap=:cividis)
Colorbar(g3[1, 2], s, ticklabelsize=tick_label_size,)
colgap!(g3, bar_gap)

ax = Axis(f[2,3], yticklabelsize=tick_label_size, xticklabelsize=tick_label_size,
          xlabel="Normal score of de-trended Tmax", xlabelsize=axis_label_size,
          title="F", titlesize=title_size, titlealign=:left)
ylims!(ax, 0, 0.61)
hist!(f[2,3], norm_stat_gr.dt_tmax, normalization=:pdf, color=wong[3])

vario = EmpiricalVariogram(norm_stat_gr, :dt_tmax, maxlag=KRIG_RANGE,
                           nlags=NLAGS,)
gamma = GeoStatsFunctions.fit(PentasphericalVariogram, vario, h->1)

dem_fn = "data/ca_dem_1km.tif"
dem = get_ras_arr(dem_fn)
ras_meta = get_meta(dem_fn)
gt = ras_meta.geotransform

nodata = -9999.
has_dem = dem .!= nodata
locs = findall(dem .!= nodata)

# get closest center of grid cell to each validation point
y_coord, x_coord = get_y_x_coords(ras_meta)
grid_locs = hcat([x_coord[x[2]] for x in findall(has_dem)],
                 [y_coord[x[1]] for x in findall(has_dem)])'
out_points = PointSet(grid_locs)
out_elev = dem[has_dem]

sol = norm_stat_gr |> InterpolateNeighbors(out_points, Kriging(gamma, 0), prob=true, neighborhood=MetricBall(KRIG_RANGE), maxneighbors=10)
arr = ones(size(dem))*NaN
arr[has_dem] = mean.(sol.dt_tmax)

g4 = f[1, 4] = GridLayout()
ax = Axis(g4[1,1], aspect=1, xticklabelsize=tick_label_size, yticklabelsize=tick_label_size,
          xticks=(-4e5:4e5:4e5, string.([-400:400:400...])),
          xlabel="Easting (km)", xlabelsize=axis_label_size, 
          yticks=(-4e5:4e5:5e5, string.([-400:400:400...])),
          limits=(-4.2e5, 5.7e5, -6.6e5, 4.9e5),
          title="G", titlesize=title_size, titlealign=:left)
ax.yticklabelsize = 0
h = heatmap!(ax, x_coord, y_coord, arr', colormap=:cividis, rasterize=2)
Colorbar(g4[1, 2], h, ticklabelsize=tick_label_size,)
colgap!(g4, bar_gap)

ax = Axis(f[2,4], yticklabelsize=tick_label_size, xticklabelsize=tick_label_size, xlabel="\"Gaussian & stationary\" estimates", xlabelsize=axis_label_size,
          title="H", titlesize=title_size, titlealign=:left)
ylims!(ax, 0, 0.61)
hist!(f[2,4], arr[has_dem], normalization=:pdf, color=wong[3])


means = georef((dt_tmax=mean.(sol.dt_tmax),), sol.geometry)
tmax_means = revert(norm_pipe, means, norm_cache)
arr[has_dem] = tmax_means.dt_tmax

g5 = f[1, 5] = GridLayout()
ax = Axis(g5[1,1], aspect=1, xticklabelsize=tick_label_size, yticklabelsize=tick_label_size,
          xticks=(-4e5:4e5:4e5, string.([-400:400:400...])),
          xlabel="Easting (km)", xlabelsize=axis_label_size, 
          yticks=(-4e5:4e5:5e5, string.([-400:400:400...])),
          limits=(-4.2e5, 5.7e5, -6.6e5, 4.9e5),
          title="I", titlesize=title_size, titlealign=:left)
ax.yticklabelsize = 0
h = heatmap!(ax, x_coord, y_coord, arr', colormap=:viridis, rasterize=2)
Colorbar(g5[1, 2], h, ticklabelsize=tick_label_size,)
colgap!(g5, bar_gap)

ax = Axis(f[2,5], yticklabelsize=tick_label_size, xticklabelsize=tick_label_size, xlabel="Estimated \"stationary\" Tmax [C]", xlabelsize=axis_label_size,
          title="J", titlesize=title_size, titlealign=:left)
ylims!(ax, 0, 0.22)
hist!(f[2,5], arr[has_dem], normalization=:pdf, color=wong[2])


trend = ols_trend_apply(grid_locs[1, :], grid_locs[2, :], out_elev,
                        df.x, df.y, df.elev, df.tmax
                        ; MAXLAG=TREND_RANGE, MIN_TREND_STATS=MIN_TREND_STATS)

arr[has_dem] = tmax_means.dt_tmax + trend

g6 = f[1, 6] = GridLayout()
ax = Axis(g6[1,1], aspect=1, xticklabelsize=tick_label_size, yticklabelsize=tick_label_size,
          xticks=(-4e5:4e5:4e5, string.([-400:400:400...])),
          xlabel="Easting (km)", xlabelsize=axis_label_size, 
          yticks=(-4e5:4e5:5e5, string.([-400:400:400...])),
          limits=(-4.2e5, 5.7e5, -6.6e5, 4.9e5),
          title="K", titlesize=title_size, titlealign=:left)
ax.yticklabelsize = 0
h = heatmap!(ax, x_coord, y_coord, arr', colormap=:inferno, rasterize=2)
Colorbar(g6[1, 2], h, ticklabelsize=tick_label_size,)
colgap!(g6, bar_gap)

ax = Axis(f[2,6], yticklabelsize=tick_label_size, xticklabelsize=tick_label_size, xlabel="Estimated Tmax [C]", xlabelsize=axis_label_size,
          title="L", titlesize=title_size, titlealign=:left)
ylims!(ax, 0, 0.082)
hist!(f[2,6], arr[has_dem], normalization=:pdf, color=wong[1])

save("figs/fig02.pdf", f)
f
