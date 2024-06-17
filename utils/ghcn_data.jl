using NCDatasets, DataFrames, ArchGDAL, Proj

const AG = ArchGDAL

include("ras_utils.jl")

function get_grid_center_elev(x, y ; grid="cimis_6km", NODATA=-9999.)
    dem_fn = "data/ca_dem_$grid.tif"
    dem = get_ras_arr(dem_fn)
    ras_meta = get_meta(dem_fn)
    gt = ras_meta.geotransform

    xind = Int.((x .- gt[1]) .÷ gt[2]) .+ 1
    yind = Int.((y .- gt[4]) .÷ gt[6]) .+ 1

    outside_dem = xind .> size(dem, 2) .|| xind .< 1 .|| yind .> size(dem, 1) .|| yind .< 1
    inside_dem = .!outside_dem

    elev = ones(length(xind)).*NODATA
    elev[inside_dem] = dem[CartesianIndex.(yind[inside_dem], xind[inside_dem])]

    xcoord = gt[1] .+ gt[2] .* xind .- gt[2]/2
    ycoord = gt[4] .+ gt[6] .* yind .- gt[6]/2

    xind, yind, xcoord, ycoord, elev
end

function check_in_bound(lon, lat, area)
    if area == :ca
        fn = "data/ca_bounds/ca_bounds_4326.shp"
    elseif area == :buff
        fn = "data/ca_buff_bounds/ca_buff.shp"
    end

    bound = AG.read(fn) do f
        AG.getlayer(f) do l
            AG.getfeature(l, 0) do f
                AG.getgeom(f, 0)
            end
        end
    end
    pts = AG.createpoint.(lon, lat)
    AG.contains.((bound,), pts)
end

function north_am(yrs ; iter=1, area=:ca, proj_coord=true, grid="cimis_6km")
    year = [yrs...][iter]
    ds = NCDataset("data/daymet_v4_stnxval_tmax_na_$year.nc")

    # id string delimited with '\0'
    ids = strip.(String.(eachcol(ds["station_id"][:, :])), '\0')
    lat = ds["stn_lat"][:]
    lon = ds["stn_lon"][:]

    good_stats = fill(true, length(ids))
    good_stats = good_stats .& check_in_bound(lon, lat, :buff)

    ids = ids[good_stats]
    ids = convert(Vector{String}, ids)

    lon = lon[good_stats]
    lat = lat[good_stats]

    # removing name for now, don't need
    #name = strip.(String.(eachcol(ds["station_name"][:, :])))[good_stats]
    elev = ds["stnz"][:][good_stats]
    tmax = ds["obs"][:, :][:, good_stats]

    df = DataFrame(id=ids, lon=lon, lat=lat, elev=elev)

    tr = Proj.Transformation("EPSG:4326", "EPSG:3310")
    p_coord = tr.(df.lat, df.lon)
    df.x = [z[1] for z in p_coord]
    df.y = [z[2] for z in p_coord]

    xind, yind, gridx, gridy, grid_elev = get_grid_center_elev(df.x, df.y ; grid=grid)
    df.xind = xind
    df.yind = yind
    df.gridx = gridx
    df.gridy = gridy
    df.grid_elev = grid_elev

    if iter == length(yrs)
        # remove non-uniqe
        unq = .!nonunique(df[:, [:x, :y]])
        tmax = tmax[:, unq]
        df = df[unq, :]
        return df, tmax
    else
        prior_df, prior_tmax = north_am(yrs ; iter=iter+1, area=area)
        current_in_prior = in.(ids, (prior_df.id,))
        df = df[current_in_prior, :]
        tmax = tmax[:, current_in_prior]

        prior_in_current = in.(prior_df.id, (ids,))
        prior_tmax = prior_tmax[:, prior_in_current]
        #df = prior_df[prior_in_current, :]

        return df, vcat( tmax, prior_tmax,)
    end
end
