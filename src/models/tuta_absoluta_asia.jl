# Guimapi, Srinivasan, Tonnang, Sotelo-Cardona & Mohamed (2020): the natural spread of the tomato leafminer
# Tuta absoluta across Asia from its first record in Turkey, with the cellular automaton of Guimapi et al.
# (2016) driven by monthly temperature and relative humidity. A land cell reached by the invasion becomes
# invaded in a month warm and humid enough for the pest, and exposed (crossed but not established)
# otherwise; an exposed cell becomes invaded in a later suitable month. Both pass the invasion on.

const TUTA_S = Int8(0)
const TUTA_E = Int8(1)
const TUTA_I = Int8(2)

"""
    build_tuta_absoluta_asia(; radius=1) -> CAModel

The cellular automaton of [Guimapi2016](@cite) as [Guimapi2020](@cite) used it
for Asia, on monthly steps: each land cell (`land`) is susceptible (0),
exposed (1) or invaded (2) (`state`). In each month, in every cell at once,

- a susceptible cell with an exposed or invaded cell among its Moore
  neighbours within `radius` cells becomes invaded if its mean temperature
  that month is at least `temperature` (°C) and its mean relative humidity at
  least `humidity` (%), and exposed otherwise;
- an exposed cell becomes invaded in a month that meets both thresholds;
- an invaded cell stays invaded.

The covariates `temperature` and `humidity` hold a slice per month. The
layers `reached` and `invaded` record the month in which a cell was first
reached and invaded (counting from 1; -1 if never). The vegetation index and
tomato production of [Guimapi2016](@cite) are left out, as in
[Guimapi2020](@cite).
"""
function build_tuta_absoluta_asia(; radius::Integer=1)
    radius >= 1 || throw(ArgumentError("radius must be at least 1, got $radius"))
    rule = LocalRule(Moore(radius)) do c, nbrs, ctx
        p = ctx.p
        (c.land && c.state != TUTA_I) || return c
        ok = covariate(ctx, :temperature) >= p.temperature &&
             covariate(ctx, :humidity) >= p.humidity
        month = Int16(ctx.t + 1)
        if c.state == TUTA_E
            return ok ? (; c..., state=TUTA_I, invaded=month) : c
        end
        any(n -> n.state != TUTA_S, nbrs) || return c
        return ok ? (; c..., state=TUTA_I, reached=month, invaded=month) :
               (; c..., state=TUTA_E, reached=month)
    end
    states = [StateVariable(:state, TUTA_S; domain=Int8(0):Int8(2),
                            description="0 susceptible, 1 exposed, 2 invaded"),
              StateVariable(:land, false; description="the cell is land in the study area"),
              StateVariable(:reached, Int16(-1);
                            description="month in which the cell was first exposed or invaded (-1 if never)"),
              StateVariable(:invaded, Int16(-1);
                            description="month in which the cell was invaded (-1 if never)")]
    parameters = [Parameter(:temperature, 22.0;
                            description="lowest monthly mean temperature for invasion (°C)"),
                  Parameter(:humidity, 55.0;
                            description="lowest monthly mean relative humidity for invasion (%)")]
    return CAModel(rule; name=:tuta_absoluta_asia, states, parameters,
                   covariates=(:temperature, :humidity),
                   description="Guimapi et al. (2020) natural spread of Tuta absoluta across Asia")
end

const _TUTA_CLIMATE = Ref{Any}(nothing)

"""
    tuta_climate() -> NamedTuple

The monthly climate of Asia bundled with the model: CRU CL 2.0 mean
temperature (°C) and relative humidity (%), 1961–1990
([New2002](@cite)), averaged over cells of 1/3 degree on a grid of 267 rows
from 78° N and 462 columns from 26° E. Returns `temperature` and `humidity`
(`267 × 462 × 12` arrays, January first, `NaN` off the land of the study
area), `land` (a `BitMatrix`) and the latitudes and longitudes of the cell
centres (`lat`, `lon`). The data are under the Open Database License (see
`models/tuta_absoluta_asia/data/NOTICE.md`).
"""
function tuta_climate()
    _TUTA_CLIMATE[] === nothing || return _TUTA_CLIMATE[]
    path = joinpath(MODELS_DIR, "tuta_absoluta_asia", "data", "cru_asia.bin")
    raw = Vector{Int16}(undef, filesize(path) ÷ 2)
    read!(path, raw)
    rec = reshape(ltoh.(raw), 26, :)
    n1, n2 = 267, 462
    temperature = fill(NaN, n1, n2, 12)
    humidity = fill(NaN, n1, n2, 12)
    land = falses(n1, n2)
    for k in axes(rec, 2)
        i, j = Int(rec[1, k]), Int(rec[2, k])
        land[i, j] = true
        temperature[i, j, :] .= rec[3:14, k] ./ 10
        humidity[i, j, :] .= rec[15:26, k] ./ 10
    end
    lat = [78 - (i - 0.5) / 3 for i in 1:n1]
    lon = [26 + (j - 0.5) / 3 for j in 1:n2]
    _TUTA_CLIMATE[] = (; temperature, humidity, land, lat, lon)
    return _TUTA_CLIMATE[]
end

"""
    tuta_cell(lat, lon; within=3) -> CartesianIndex{2}

The cell of [`tuta_climate`](@ref)'s grid holding the point (degrees north and
east), or the nearest land cell within `within` cells of it if that cell is
not land.
"""
function tuta_cell(lat::Real, lon::Real; within::Integer=3)
    c = tuta_climate()
    i = floor(Int, (78 - lat) * 3) + 1
    j = floor(Int, (lon - 26) * 3) + 1
    (1 <= i <= size(c.land, 1) && 1 <= j <= size(c.land, 2)) ||
        throw(ArgumentError("($lat, $lon) lies outside the grid (11° S–78° N, 26° E–180° E)"))
    c.land[i, j] && return CartesianIndex(i, j)
    best = nothing
    for a in max(1, i - within):min(size(c.land, 1), i + within),
        b in max(1, j - within):min(size(c.land, 2), j + within)
        c.land[a, b] || continue
        d = hypot(a - i, b - j)
        (best === nothing || d < best[1]) && (best = (d, CartesianIndex(a, b)))
    end
    best === nothing &&
        throw(ArgumentError("no land of the study area within $within cells of ($lat, $lon)"))
    return best[2]
end

"""
    problem_tuta_absoluta_asia(; radius=1, origin=(38.4, 27.1), start_month=1, years=25,
                               seed=0, pmap=(;), kwargs...) -> CAProblem

The setup of [Guimapi2020](@cite): Asia on the grid of [`tuta_climate`](@ref),
susceptible everywhere but the cell at `origin` (degrees north and east; by
default İzmir, near the first record in Turkey, in 2009), which is invaded,
for `years` years of monthly steps from the month `start_month`. The
covariates repeat the climatology every 12 steps. The model is deterministic.
"""
function problem_tuta_absoluta_asia(; radius::Integer=1, origin=(38.4, 27.1),
                                    start_month::Integer=1, years::Integer=25,
                                    seed::Integer=0, pmap=(;), kwargs...)
    1 <= start_month <= 12 ||
        throw(ArgumentError("start_month must lie in 1:12, got $start_month"))
    years >= 1 || throw(ArgumentError("years must be at least 1, got $years"))
    c = tuta_climate()
    n1, n2 = size(c.land)
    order = [mod1(start_month + k - 1, 12) for k in 1:12]
    temperature = c.temperature[:, :, order]
    humidity = c.humidity[:, :, order]
    O = tuta_cell(origin...)
    state = fill(TUTA_S, n1, n2)
    state[O] = TUTA_I
    reached = fill(Int16(-1), n1, n2)
    invaded = fill(Int16(-1), n1, n2)
    reached[O] = invaded[O] = 0
    model = build_tuta_absoluta_asia(; radius)
    return CAProblem(model, (; state, land=copy(c.land), reached, invaded), 12 * years,
                     pmap;
                     space=Grid(n1, n2), seed, covariates=(; temperature, humidity),
                     kwargs...)
end

"""
    tuta_arrival(u, lat, lon) -> NamedTuple

The years after the start at which the invasion first reached and invaded
the cell at `(lat, lon)` (see [`tuta_cell`](@ref)) in the state `u`, or
`nothing` for either that has not happened.
"""
function tuta_arrival(u, lat::Real, lon::Real)
    I = tuta_cell(lat, lon)
    years(m) = m < 0 ? nothing : m / 12
    return (; reached=years(u.reached[I]), invaded=years(u.invaded[I]))
end

_register!(:tuta_absoluta_asia; build=build_tuta_absoluta_asia,
           problem=problem_tuta_absoluta_asia)
