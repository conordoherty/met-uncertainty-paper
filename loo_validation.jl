using GeoStats
using JLD2
using ProgressMeter
using Random

Random.seed!(1)

include("utils/ghcn_data.jl")
include("utils/ols_trend.jl")
include("utils/result.jl")

TREND_RANGE = 1.00e5
KRIG_RANGE = 1.00e5
GRID = "1km"

MIN_TREND_STATS = 10
NLAGS = 12
MAX_NEIGH = 100
NODATA = -9999.

quants = [0.005:0.005:0.995...]
ints = [.99:-.01:.01...]

df, tmax = north_am(2022 ; grid=GRID)

function run_location(ind, doy)
    # tmax for all stations on DOY
    data = tmax[doy, :]

    # skip if no data 
    ismissing(data[ind]) && return NODATA

    # skip if not in dem
    df[ind, :grid_elev] == NODATA && return NODATA

    # get row of prediction location
    pred_stat = df[ind:ind, :]
    true_tmax = data[ind]

    # drop stats with no data
    has_data = (!ismissing).(data)
    df_stats = df[has_data, :]
    data = data[has_data]

    # remove prediction location
    new_ind = findall(df_stats.id .== pred_stat.id)[1]
    df_stats = df_stats[Not(new_ind), :]
    data = data[Not(new_ind)]

    # skip if prediction location in the same grid cell as another
    # station with data
    same_coord = innerjoin(pred_stat[:, [:gridx, :gridy]],
                           df_stats[:, [:gridx, :gridy]],
                           on=[:gridx, :gridy])
    size(same_coord, 1) >= 1 && return NODATA

    # detrend tmax, store in the df
    df_stats.tmax = data
    df_stats.dt_tmax, _ = ols_trend_domain(df_stats.x, df_stats.y,
                                           df_stats.elev, data
                                           ; MAXLAG=TREND_RANGE,
                                           MIN_TREND_STATS=MIN_TREND_STATS)

    # normal score transform
    stat_gr = georef(df_stats[:, [:x, :y, :dt_tmax]], (:x, :y))
    norm_pipe = Quantile()
    norm_stat_gr, norm_cache = apply(norm_pipe, stat_gr)

    vario = EmpiricalVariogram(norm_stat_gr, :dt_tmax,
                               maxlag=KRIG_RANGE, nlags=NLAGS,)
    gamma = GeoStatsFunctions.fit(PentasphericalVariogram, vario, h->1)
 
    out_points = PointSet(hcat(pred_stat.gridx, pred_stat.gridy)')

    sol = norm_stat_gr |>
          InterpolateNeighbors(out_points, Kriging(gamma), prob=true,
                               neighborhood=MetricBall(KRIG_RANGE),
                               maxneighbors=MAX_NEIGH, minneighbors=2)

    # if can't produce estimate at location (eg # neighbors 
    # < minneighbors) it returns missing
    # will error when trying to revert normal score transform
    if ismissing(sol.dt_tmax[1])
        println(ind, " ", doy)
    end

    # normal pred distribution quantiles
    norm_quants = quantile.(sol.dt_tmax, (quants,))
    norm_quants = reduce(hcat, norm_quants)'

    function denorm_col(ind)
        gr = georef((dt_tmax=norm_quants[:, ind],), sol.geometry)
        revert(norm_pipe, gr, norm_cache).dt_tmax
    end
    
    denorm_quants = reduce(hcat, [denorm_col(i) for i in 1:length(quants)])

    # calculate trend at prediction location
    trend = ols_trend_apply(pred_stat.gridx, pred_stat.gridy,
                            pred_stat.grid_elev,
                            df_stats.x, df_stats.y, df_stats.elev,
                            df_stats.tmax ; MAXLAG=TREND_RANGE,
                            MIN_TREND_STATS=MIN_TREND_STATS)
    pred_stat.trend = trend

    # add trend back to pred distribution quantiles
    detrend_quants = denorm_quants .+ trend
 
    Result(vario, gamma, df_stats, data, true_tmax, pred_stat, sol,
           detrend_quants[1, :])
end

@showprogress Threads.@threads for DOY = 1:365
    res = [run_location(i, DOY) for i = 1:size(df, 1)]
    jldsave("output/res_$DOY.jld2", res=res)
end
