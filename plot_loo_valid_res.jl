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

median_doy = findfirst(mae_loc .== median(mae_loc))
println("median MAE DOY: ", median_doy)
println("median MAE: ", mae_loc[median_doy])
worst_pred_int_median_ind = argmax(abs.(int_day_loc[:, median_doy]-ints))
println("median worst pred int: ", ints[worst_pred_int_median_ind])
println("median worst pred: ", int_day_loc[worst_pred_int_median_ind, median_doy])
println("")

low_doy = argmin(bias_loc)
println("low worst DOY: ", low_doy)
println("low worst bias: ", bias_loc[low_doy])
worst_pred_int_low_ind = argmin(int_day_loc[:, low_doy]-ints)
println("low worst pred int: ", ints[worst_pred_int_low_ind])
println("low worst pred: ", int_day_loc[worst_pred_int_low_ind, low_doy])
println("")

high_doy = argmax(bias_loc)
println("high worst DOY: ", high_doy)
println("high worst bias: ", bias_loc[high_doy])
worst_pred_int_high_ind = argmin(int_day_loc[:, high_doy]-ints)
println("high worst pred int: ", ints[worst_pred_int_high_ind])
println("high worst pred: ", int_day_loc[worst_pred_int_high_ind, high_doy])
