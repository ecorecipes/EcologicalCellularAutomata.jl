# Builds cru_asia.bin from the CRU CL 2.0 10' climatology (New et al. 2002), files grid_10min_tmp.dat.gz
# and grid_10min_reh.dat.gz from https://crudata.uea.ac.uk/cru/data/hrg/tmc/ (Open Database License).
# Not run by the package; kept to document how the bundled file was made.
#
#   julia make_cru_asia.jl grid_10min_tmp.dat.gz grid_10min_reh.dat.gz cru_asia.bin
#
# Cells of 1/3 degree over 26-180 E, 11 S-78 N (267 rows from the north, 462 columns from the west), each
# the mean of its four 10' cells that are land, kept if at least two are land and the cell lies in Asia
# (the rough boundary in `asia`). One record per kept cell, 26 little-endian Int16: row, column, the 12
# monthly mean temperatures (0.1 degC) and the 12 monthly mean relative humidities (0.1 %), January first.
using Statistics

const LON0, LAT1 = 26.0, 78.0
const ROWS, COLS = 267, 462

function readgrid(path)
    out = Dict{Tuple{Int,Int},Vector{Float64}}()
    for l in eachline(pipeline(`gzip -dc $path`))
        f = parse.(Float64, split(l))
        lat, lon = f[1], f[2]
        i = floor(Int, (LAT1 - lat) * 6) + 1
        j = floor(Int, (lon - LON0) * 6) + 1
        (1 <= i <= 2ROWS && 1 <= j <= 2COLS) || continue
        out[(i, j)] = f[3:14]
    end
    return out
end

# A rough Asia: no Europe (west of the Turkish straits, north of the Black Sea and the Caucasus, west of
# the Ural river and mountains), no Africa (west of Suez and of the Red Sea, and the Horn) and no New Guinea.
function asia(lat, lon)
    lon < 32.5 && lat < 37 && return false
    # Africa west of a line along the Red Sea from Suez (30 N, 32.3 E) to Bab-el-Mandeb (12.6 N, 43.4 E),
    # and the Horn of Africa south of the Gulf of Aden
    lat <= 30 && lon < 32.3 + (30 - lat) * (43.4 - 32.3) / (30 - 12.6) && return false
    lat < 12.6 && lon < 52.0 && return false
    lat > 40.9 && lon < 29.5 && return false
    lat > 41.2 && lon < 41.5 && return false
    lat > 44.0 && lon < 47.5 && return false
    lat > 47.0 && lon < 52.0 && return false
    lat > 51.5 && lon < 60.0 && return false
    lon > 140.9 && lat < 0 && return false
    return true
end

function main(tmpfile, rehfile, out)
    T = readgrid(tmpfile)
    H = readgrid(rehfile)
    records = Int16[]
    for i in 1:ROWS, j in 1:COLS
        sub = [(a, b) for a in (2i - 1):(2i) for b in (2j - 1):(2j)]
        land = [s for s in sub if haskey(T, s) && haskey(H, s)]
        length(land) >= 2 || continue
        asia(LAT1 - (i - 0.5) / 3, LON0 + (j - 0.5) / 3) || continue
        t = round.(Int16, 10 .* mean(T[s] for s in land))
        h = round.(Int16, 10 .* mean(H[s] for s in land))
        append!(records, Int16[i, j], t, h)
    end
    open(out, "w") do io
        return write(io, htol.(records))
    end
    return println("wrote ", length(records) ÷ 26, " cells to ", out)
end

main(ARGS...)
