using CairoMakie, GeoIO

include("utils/ghcn_data.jl")

bound_fn = "data/ca_bounds/ca_bounds_4326.shp"
bound = GeoIO.load(bound_fn)

df, tmax = north_am(2022 ; grid="1km")

f = Figure(size=(1100,500))
ax = Axis(f[1,1], aspect=1, ylabel="Latitude",
          xlabel="Longitude", title="A", titlealign=:left)

geoms = bound.geometry[1].geoms

g(geom) = map(x->(x.lon.val, x.lat.val),
              getfield.(geom.rings[1].vertices, :coords))

for geom in geoms
    lines!(ax, g(geom), color=:black, label="Study area")
end

scatter!(ax, df.lon, df.lat, markersize=5,
         label="Weather stations")
axislegend(ax, merge=true)

gr = f[1,2] = GridLayout()
colsize!(gr, 1, 530)
ax = Axis(gr[1,1], ylabel="Tmax [°C]", xlabel="Day of year",
         title="B", titlealign=:left)
lines!(ax, tmax[:, 300], label="Coastal")
lines!(ax, tmax[:, 386], label="Inland")
lines!(ax, tmax[:, 572], label="Mountain")
axislegend(ax,)

save("figs/fig01.pdf", f)
#f
