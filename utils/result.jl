using GeoStats, DataFrames

# location result
struct Result
    vario::EmpiricalVariogram
    gamma::Variogram
    data_locs::DataFrame
    data::Vector{Float64}
    true_val::Float64
    pred_stat::DataFrame
    pred_gauss::GeoTable
    quants::Vector{Float64}
end
