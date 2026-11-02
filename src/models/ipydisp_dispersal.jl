# Agboka et al. (2024): iPyDisp, a cellular automaton for the dispersal of biological control agents.
# Cells are the centroids of a longitude-latitude grid over a region; a cell's neighbours are the cells
# within a travel distance (great-circle). Each step, every cell reached by an exposed or colonised
# neighbour becomes colonised if it meets every environmental constraint, and exposed otherwise; exposed
# cells spread too. The model is deterministic, so it reproduces the software's runs exactly.

"""
    build_ipydisp_dispersal() -> CAModel

The dispersal automaton of iPyDisp [Agboka2024](@cite) on a graph of cells
(a `GraphSpace` from [`ipydisp_neighbours`](@ref)). Each cell holds its
`status` (0 unexposed, 1 exposed, 2 colonised) and whether it is `suitable`
(it meets every constraint, see [`ipydisp_suitable`](@ref)). In each step,
in every cell at once, a cell that is not colonised and has an exposed or
colonised neighbour becomes colonised if it is suitable and exposed if not;
other cells keep their status.
"""
function build_ipydisp_dispersal()
    rule = LocalRule(GraphNeighbours()) do c, nbrs, ctx
        c.status == 2 && return c
        reached = false
        for n in nbrs
            n.status >= 1 && (reached=true; break)
        end
        reached || return c
        return (; status=c.suitable ? Int8(2) : Int8(1), suitable=c.suitable)
    end
    states = [StateVariable(:status, Int8(0);
                            description="0 unexposed, 1 exposed, 2 colonised"),
              StateVariable(:suitable, false;
                            description="the cell meets every environmental constraint")]
    return CAModel(rule; name=:ipydisp_dispersal, states,
                   description="Agboka et al. (2024) iPyDisp dispersal of biological control agents")
end

"""
    ipydisp_distance(lat1, lon1, lat2, lon2) -> Float64

The great-circle (haversine) distance in km between two points given in
degrees, on a sphere of radius 6371 km, as iPyDisp computes it.
"""
function ipydisp_distance(lat1, lon1, lat2, lon2)
    dlat = deg2rad(lat2 - lat1)
    dlon = deg2rad(lon2 - lon1)
    a = sin(dlat / 2)^2 + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * sin(dlon / 2)^2
    return 6371 * 2 * atan(sqrt(a), sqrt(1 - a))
end

"""
    ipydisp_grid(; west, east, south, north, cellsize=1, inside=(lat, lon) -> true)
        -> (; lat, lon, row, col)

The cells of iPyDisp [Agboka2024](@cite): a grid over the bounding box with
cells of `cellsize × 0.00833` degrees (about `cellsize` km), numbered row by
row from the north-west corner, cells at the east and south edges cut at the
box. A cell is kept, at the centroid of its (cut) square, if `inside` holds
there; iPyDisp keeps the cells meeting the region's polygon and places them at
the centroid of their part inside it. Returns the kept cells' coordinates and
their row and column in the grid.
"""
function ipydisp_grid(; west::Real, east::Real, south::Real, north::Real,
                      cellsize::Real=1, inside=(lat, lon) -> true)
    west < east && south < north ||
        throw(ArgumentError("the box must have west < east and south < north"))
    cellsize > 0 || throw(ArgumentError("cellsize must be positive, got $cellsize"))
    d = cellsize * 0.00833
    nx = ceil(Int, (east - west) / d)
    ny = ceil(Int, (north - south) / d)
    lat = Float64[]
    lon = Float64[]
    row = Int[]
    col = Int[]
    for i in 0:(ny - 1), j in 0:(nx - 1)
        x0, x1 = min(west + d * j, east), min(west + d * (j + 1), east)
        y0, y1 = max(north - d * i, south), max(north - d * (i + 1), south)
        y, x = (y0 + y1) / 2, (x0 + x1) / 2
        inside(y, x) || continue
        push!(lat, y)
        push!(lon, x)
        push!(row, i + 1)
        push!(col, j + 1)
    end
    return (; lat, lon, row, col)
end

"""
    ipydisp_neighbours(lat, lon, travel) -> Vector{Vector{Int}}

The neighbourhood of iPyDisp [Agboka2024](@cite): for each cell, the other
cells whose centroids lie less than `travel` km away (see
[`ipydisp_distance`](@ref)), in increasing order. Pass it to `GraphSpace`.
"""
function ipydisp_neighbours(lat::AbstractVector, lon::AbstractVector, travel::Real)
    length(lat) == length(lon) ||
        throw(ArgumentError("lat and lon must have the same length"))
    n = length(lat)
    return [[k
             for k in 1:n
             if k != i && ipydisp_distance(lat[i], lon[i], lat[k], lon[k]) < travel]
            for i in 1:n]
end

"""
    ipydisp_suitable(constraints, n) -> Vector{Bool}

Whether each of `n` cells meets every constraint of iPyDisp
[Agboka2024](@cite). `constraints` is a collection of `(values, lower, upper)`,
one per environmental layer, with a value per cell; a cell meets a constraint
when `lower < value <= upper`. With no constraints every cell is suitable (the
software's default bounds, -1 and 1, around a layer of zeros).
"""
function ipydisp_suitable(constraints, n::Integer)
    ok = fill(true, n)
    for (values, lower, upper) in constraints
        length(values) == n ||
            throw(ArgumentError("a constraint has $(length(values)) values for $n cells"))
        ok .&= (lower .< values) .& (values .<= upper)
    end
    return ok
end

"""
    ipydisp_example(; cellsize=10, travel=15, seed=0) -> NamedTuple

A synthetic setting for iPyDisp: a box from 34° to 36° E and 1° S to 1° N
(roughly central Kenya) cut to an ellipse, with one constraint, a smooth
random layer in [0, 1] drawn from `seed` (a random cluster landscape,
averaged over each cell's Moore neighbourhood), that cells must exceed 0.3 to
be colonised; travel of 15 km joins each cell of about 9.3 km to its eight
neighbours. Returns the `grid`, the `constraints` and the `neighbours`.
"""
function ipydisp_example(; cellsize::Real=10, travel::Real=15, seed::Integer=0)
    inside(lat, lon) = (lon - 35)^2 + lat^2 <= 1
    g = ipydisp_grid(; west=34, east=36, south=-1, north=1, cellsize, inside)
    rng = Random.Xoshiro(seed)
    m = maximum(g.row)
    n = maximum(g.col)
    L = random_cluster_landscape(rng, max(m, n), 0.4)
    S = [begin
             is = max(1, i - 1):min(size(L, 1), i + 1)
             js = max(1, j - 1):min(size(L, 2), j + 1)
             sum(L[is, js]) / (length(is) * length(js))
         end
         for i in 1:m, j in 1:n]
    S = (S .- minimum(S)) ./ (maximum(S) - minimum(S))
    values = [S[g.row[k], g.col[k]] for k in eachindex(g.lat)]
    return (; grid=g, constraints=[(values, 0.3, 1.0)],
            neighbours=ipydisp_neighbours(g.lat, g.lon, travel))
end

"""
    problem_ipydisp_dispersal(; grid=nothing, neighbours=nothing, constraints=nothing,
                              start=nothing, steps=24, seed=0, kwargs...) -> CAProblem

A run of iPyDisp [Agboka2024](@cite) for `steps` time steps (months or years,
as the user reads them): the cells in `start` (indices into the grid; iPyDisp
takes the cells nearest the released sites) begin colonised. Without a grid,
the setting is [`ipydisp_example`](@ref) drawn from `seed`; without `start`,
the run starts from the cell nearest the grid's mean position. The model is
deterministic.
"""
function problem_ipydisp_dispersal(; grid=nothing, neighbours=nothing, constraints=nothing,
                                   start=nothing, steps::Integer=24, seed::Integer=0,
                                   kwargs...)
    if grid === nothing
        ex = ipydisp_example(; seed)
        grid, neighbours = ex.grid, ex.neighbours
        constraints === nothing && (constraints = ex.constraints)
    end
    neighbours === nothing &&
        throw(ArgumentError("give the neighbourhood (ipydisp_neighbours) with the grid"))
    n = length(grid.lat)
    length(neighbours) == n ||
        throw(ArgumentError("the neighbourhood has $(length(neighbours)) cells, the grid $n"))
    suitable = constraints === nothing || isempty(constraints) ? fill(true, n) :
               ipydisp_suitable(constraints, n)
    if start === nothing
        clat, clon = sum(grid.lat) / n, sum(grid.lon) / n
        start = [argmin([ipydisp_distance(grid.lat[k], grid.lon[k], clat, clon)
                         for k in 1:n])]
    end
    all(k -> 1 <= k <= n, start) || throw(ArgumentError("start cells must lie in 1:$n"))
    status = zeros(Int8, n)
    status[start] .= 2
    return CAProblem(build_ipydisp_dispersal(), (; status, suitable), steps, (;);
                     space=GraphSpace(neighbours), seed, kwargs...)
end

"""
    ipydisp_counts(u) -> NamedTuple

The numbers of `unexposed`, `exposed` and `colonised` cells in a state.
"""
function ipydisp_counts(u)
    return (; unexposed=count(==(0), u.status), exposed=count(==(1), u.status),
            colonised=count(==(2), u.status))
end

_register!(:ipydisp_dispersal; build=build_ipydisp_dispersal,
           problem=problem_ipydisp_dispersal)
