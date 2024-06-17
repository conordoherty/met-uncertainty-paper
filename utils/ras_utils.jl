using ArchGDAL

const AG = ArchGDAL

struct RasterMeta
    geotransform::NTuple{6, Float64}
    height::Int
    width::Int
    EPSG_code::Int
end

struct Raster
    fn::String
    height::Int
    width::Int
    nodata::Number
    data::Vector{Number}
    data_locs::Vector{CartesianIndex{2}}
    gt::Vector{AbstractFloat}
    xcoord::Vector{AbstractFloat}
    ycoord::Vector{AbstractFloat}
end

function Raster(fn)
    ras = AG.readraster(fn)
    height = AG.height(ras)
    width = AG.width(ras)
    nodata = AG.getnodatavalue(AG.getband(ras, 1))
    arr = ras[:, :, 1]

    # need to set length of cart inds for some reason
    data_locs = findall(arr .!= nodata)

    data = arr[data_locs]
    gt = AG.getgeotransform(ras)

    # offset by half grid cell relative to gt
    xind = [range(start=gt[1], step=gt[2], length=width)...] .+ gt[2]/2
    yind = [range(start=gt[4], step=gt[6], length=height)...] .+ gt[6]/2

    Raster(fn, height, width, nodata, data, data_locs, gt, xind, yind)
end

function plot_arr(ras::Raster, data::Vector)
    arr = ones(ras.width, ras.height) .* NaN
    arr[ras.data_locs] = data
    arr
end

function write_raster(data_arr::Matrix{Float64}, out_filename::String, ras_meta::RasterMeta, nodata=-9999.)
    AG.create(
        out_filename,
        driver=AG.getdriver("GTiff"),
        width=ras_meta.width,
        height=ras_meta.height,
        nbands=1,
        dtype=Float64
    ) do dataset
        AG.write!(dataset, data_arr, 1)
        AG.setgeotransform!(dataset, [ras_meta.geotransform...])
        AG.setproj!(dataset, AG.toWKT(AG.importEPSG(ras_meta.EPSG_code)))
        AG.getband(dataset, 1) do rasterband
            AG.setnodatavalue!(rasterband, nodata)
        end
    end
end

function get_meta(fn::String)
    ras_meta = AG.read(fn) do ds
        gt = AG.getgeotransform(ds)
        epsg = parse(Int, AG.getattrvalue(AG.importWKT(AG.getproj(ds)), "AUTHORITY", 1))
        height, width = AG.getband(ds, 1) do band
            return AG.height(band), AG.width(band)
        end
        return RasterMeta(Tuple(gt), height, width, epsg)
    end
    return ras_meta
end

function get_y_x_coords(geotransform::NTuple{6, Float64}, height::Int, width::Int)
    gt = geotransform
    y_cell_size = gt[6]
    x_cell_size = gt[2]
    ul_cell_center = (gt[4]+y_cell_size/2, gt[1]+x_cell_size/2)
    y_coords = [ul_cell_center[1]:y_cell_size:ul_cell_center[1]+(height-1)*y_cell_size...]
    x_coords = [ul_cell_center[2]:x_cell_size:ul_cell_center[2]+(width-1)*x_cell_size...]

    return y_coords, x_coords
end

get_y_x_coords(ras_meta::RasterMeta) = get_y_x_coords(ras_meta.geotransform, ras_meta.height, ras_meta.width)

function get_ras_arr(fn::String)
    # archgdal reads array as buffer (cols, rows)
    # so need to transpose
    AG.readraster(fn) do ds
       return Matrix(ds[:, :, 1]')
    end
end

# doens't do bounds checking (raster size not passed)
function get_index_ras(x_coord, y_coord, gt)
    x_ind = Integer((x_coord - gt[1]) ÷ gt[2]) + 1
    y_ind = Integer((y_coord - gt[4]) ÷ gt[6]) + 1
    CartesianIndex(y_ind, x_ind)
end

# findfirst(locs .== ind) is trying to broadcast over ind
# can't figure ^ out so wrote a stupid func instead
function get_index_locs(locs, ind)
    for i in 1:size(locs)[1]
        if locs[i] == ind
            return i
        end
    end
    # didnt' find a match
    return -1
end


#function ind_lookup

# WIP
#function downsample_by_factor(geotransform::NTuple{6, Float64}, arr::Matrix{Float64}, factor::Int ;
#                              min_pts::Int=factor^2÷2, agg_fn::Function=mean, nodata::Float64=-9999.)
#    ulx, uly = geotransform[1], geotransform[4]
#    new_dims = div.(size(arr), factor, RoundUp)
#
#    # add nodata buffer to make input dims multiple of `factor`
#    height_remainder = size(arr)[1] % factor
#    width_remainder = size(arr)[2] % factor
#    buffered_arr = hcat(arr, ones(size(arr)[1], width_remainder).*nodata)
#    buffered_arr = vcat(buffered_arr, ones(height_remainder, size(buffered_arr)[2]).*nodata)
#end
