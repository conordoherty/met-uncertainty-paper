using JLD2, Statistics
#using GLMakie
using CairoMakie

ints = [.99:-.01:.01...]
mae(x, y) = x - y .|> abs |> mean

int_day_loc = load("output/loo_valid_res.jld2")["int_day"]

bias_loc = mean(int_day_loc .- ints, dims=1)[1, :]
mae_loc = mean(abs.(int_day_loc .- ints), dims=1)[1, :]

f = Figure(size=(800, 400))

ax = Axis(f[1,1], ylabel="Prediction interval MAE", title="A", titlealign=:left)
barplot!(ax, 1:365, mae_loc ; gap=0)
xlims!(ax, -1,367)
ylims!(ax, 0, 0.043)

ax = Axis(f[2,1], ylabel="Prediction interval bias", xlabel="Day of year", title="B", titlealign=:left)
barplot!(ax, 1:365, bias_loc ; gap=0)
xlims!(ax, -1,367)

save("figs/fig04.pdf", f)

f = Figure(size=(900, 300))

ax = Axis(f[1,1], ylabel="Actual fraction within interval",
          xlabel="Theoretical prediction interval",
          title="A", titlealign=:left, aspect=1, limits=(0,1,0,1))
med_ind = findfirst(mae_loc .== median(mae_loc))
lines!(ax, ints, ints, color=:black, linestyle=:dash) 
lines!(ax, ints, int_day_loc[:, med_ind]) 

ax = Axis(f[1,2], ylabel="Actual fraction within interval", xlabel="Theoretical prediction interval",
          title="B", titlealign=:left, aspect=1, limits=(0,1,0,1))
lines!(ax, ints, ints, color=:black, linestyle=:dash) 
lines!(ax, ints, int_day_loc[:, argmax(bias_loc)]) 

ax = Axis(f[1,3], ylabel="Actual fraction within interval", xlabel="Theoretical prediction interval",
          title="C", titlealign=:left, aspect=1, limits=(0,1,0,1))
lines!(ax, ints, ints, color=:black, linestyle=:dash) 
lines!(ax, ints, int_day_loc[:, argmin(bias_loc)]) 

save("figs/fig05.pdf", f)
#f
