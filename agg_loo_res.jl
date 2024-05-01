using JLD2, DataFrames, GeoStats
using ProgressMeter

include("utils/result.jl")

NODATA = -9999.
quants = [0.005:0.005:0.995...]
ints = [.99:-.01:.01...]

int_day = zeros(99, 365)
@showprogress Threads.@threads for doy = 1:365
    res = load("output/res_$doy.jld2")["res"]
    res = res[res .!= NODATA]

	within = zeros(99)
	for i in 1:99
    	above_bottom = [res[j].quants[i] for j in 1:length(res)] .<
                       [res[j].true_val for j in 1:length(res)]
    	below_top = [res[j].quants[end-(i-1)] for j in 1:length(res)] .> 
                    [res[j].true_val for j in 1:length(res)]

    	within[i] = mean(above_bottom .& below_top)
        int_day[:, doy] = within
	end
end
jldsave("output/loo_valid_res.jld2", int_day=int_day)
