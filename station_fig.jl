using CairoMakie, GeoIO

include("utils/ghcn_data.jl")

bound_fn = "data/ca_bounds/ca_bounds_4326.shp"
bound = GeoIO.load(bound_fn)

df, _ = north_am(2022 ; grid="1km")

f = Figure(size=(600,600))
ax = Axis(f[1,1], aspect=1, ylabel="Latitude", xlabel="Longitude")

geoms = bound.geometry[1].geoms

for geom in geoms
    lines!(ax, map(x -> (x.coords[1], x.coords[2]), geom.rings[1].vertices),
           color=:black, label="Study area")
end

scatter!(ax, df.lon, df.lat, markersize=5, label="Weather stations")
axislegend(ax, merge=true)

save("figs/fig01.pdf", f)
#f
