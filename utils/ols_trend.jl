using Distances: euclidean
using LinearAlgebra: I, pinv, norm, dot
using NearestNeighbors

make_X_mat(x, y, elev) = hcat(ones(length(x)), x, y, elev,)

function ols_trend_location(x, y, elev, obs ; x_mat_fn=make_X_mat)
    X = x_mat_fn(x, y, elev)
    X \ obs
end

function ols_trend_domain(x, y, elev, obs ; MAXLAG=3.5e5, MIN_TREND_STATS=1,
                          x_mat_fn=make_X_mat)
    stat_locs = hcat(x, y)'
    tree = KDTree(stat_locs)

    dt_obs = zeros(length(obs))
    trend = zeros(length(obs))
    stat_counts = zeros(Int, length(obs))
    for j in 1:length(x)
        num_stats = max(inrangecount(tree, [x[j], y[j]], MAXLAG),
                        MIN_TREND_STATS)
        stats, _ = knn(tree, [x[j], y[j]], num_stats, true)

        x_stat = x[stats]
        y_stat = y[stats]
        elev_stat = elev[stats]
        obs_stat = obs[stats]

        X = x_mat_fn(x_stat, y_stat, elev_stat)
        β̂ = X \ obs_stat
        #trend[j] = dot([1, x[j], y[j], elev[j]], β̂)
        trend[j] = dot(X[1, :], β̂)
        dt_obs[j] = obs[j] - trend[j]
        stat_counts[j] = num_stats
    end
    dt_obs, trend, stat_counts
end

function ols_trend_apply(x_pred, y_pred, elev_pred, x_obs, y_obs, elev_obs, obs
                         ; MAXLAG=3.5e5, MIN_TREND_STATS=1, return_coefs=false)
    stat_locs = hcat(x_obs, y_obs)'
    tree = KDTree(stat_locs)

    coef_dict = Dict{Vector, Vector}()
    trend = zeros(length(x_pred))

    if return_coefs
        coef_arr = zeros(4, length(x_pred))
    end

    for i in 1:length(x_pred)
        num_stats = max(inrangecount(tree, [x_pred[i], y_pred[i]], MAXLAG),
                        MIN_TREND_STATS)
        stats, _ = knn(tree, [x_pred[i], y_pred[i]], num_stats, true)

        x_stat = x_obs[stats]
        y_stat = y_obs[stats]
        elev_stat = elev_obs[stats]
        obs_stat = obs[stats]

        if !haskey(coef_dict, stats)
            coefs = ols_trend_location(x_stat, y_stat, elev_stat, obs_stat)
            coef_dict[stats] = coefs
        else
            coefs = coef_dict[stats]
        end

        if return_coefs
            coef_arr[:, i] = coefs
        end

        loc_covariates = vcat([1], [x_pred[i], y_pred[i], elev_pred[i]])
        trend[i] = sum(loc_covariates .* coefs)
    end

    if return_coefs
        return trend, coef_arr'
    else
        return trend
    end
end
