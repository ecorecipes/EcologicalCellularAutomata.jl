# Mercader, Siegert, Liebhold & McCullough (2011): the local spread of the emerald ash borer, Agrilus
# planipennis, a coupled map lattice of 50 m cells holding ash phloem and beetles, in yearly steps. Adults
# emerge from larvae of one or two years, disperse by a negative exponential kernel, lay eggs as far as
# the phloem allows, and the larvae eat phloem and develop in one or two years depending on their density.
# Growth, feeding and development are the model's rule; dispersal, normalised source by source over a
# window and biased away from cells without phloem, is a callback stored with the problem, which leaves in
# each cell the dispersed adults that will lay there the next year.

"""
    build_emerald_ash_borer() -> CAModel

The emerald ash borer model of [Mercader2011](@cite), in yearly steps. Each
cell holds ash `phloem` (m²), the dispersed `adults` (females) that will lay
there this year, the larvae laid last year that emerge next year
(`larvae1`), the larvae that take two years, laid in odd and in even years
(`larvae2_odd`, `larvae2_even`), and whether it lies in the border around the
landscape (`border`), where dispersing beetles are lost, and its phloem at the
start (`phloem0`). Each year, in every
cell at once:

1. the adults lay `min(r0 adults, phloem / C)` larvae, which eat `C` m² of
   phloem each;
2. a share `p2 = -0.194 ln(x) + 0.986` (kept within 0 and 1) of the larvae take
   two years to develop, where `x` is the larvae per m² of the phloem left;
   the rest emerge next year.

The emergence and dispersal of adults are [`eab_disperse!`](@ref), a callback
stored with the problem.
"""
function build_emerald_ash_borer()
    rule = LocalRule(Offsets((1, 0))) do c, nbrs, ctx
        p = ctx.p
        odd = isodd(ctx.t)   # the source's iteration ctx.t + 2
        new = max(min(p.r0 * c.adults, c.phloem / p.C), 0.0)
        phloem = max(c.phloem - p.C * new, 0.0)
        p2 = clamp(-0.194 * log((new + 0.001) / (phloem + 0.0001)) + 0.986, 0.0, 1.0)
        return (; phloem, adults=0.0, larvae1=(1 - p2) * new,
                larvae2_odd=odd ? p2 * new : c.larvae2_odd,
                larvae2_even=odd ? c.larvae2_even : p2 * new, border=c.border,
                phloem0=c.phloem0)
    end
    states = [StateVariable(:phloem, 0.0; description="ash phloem (m²)"),
              StateVariable(:adults, 0.0;
                            description="dispersed adults that lay in the cell this year"),
              StateVariable(:larvae1, 0.0; description="larvae emerging next year"),
              StateVariable(:larvae2_odd, 0.0;
                            description="two-year larvae laid in odd years"),
              StateVariable(:larvae2_even, 0.0;
                            description="two-year larvae laid in even years"),
              StateVariable(:border, false;
                            description="a border cell, where dispersing beetles are lost"),
              StateVariable(:phloem0, 0.0;
                            description="the cell's phloem at the start (m²)")]
    parameters = [Parameter(:r0, 10.35; description="adult female offspring per female"),
                  Parameter(:C, 0.0113; description="phloem eaten by a larva (m²)"),
                  Parameter(:b, 0.037; description="dispersal rate of the kernel (1/m)"),
                  Parameter(:cell, 50.0; description="side of a cell (m)"),
                  Parameter(:window, 50;
                            description="half-width of the dispersal window (cells)"),
                  Parameter(:self_distance, 0.25;
                            description="distance of a cell from itself (cells)"),
                  Parameter(:behaviour, :avoidance;
                            description=":standard, :avoidance or :resource_bias"),
                  Parameter(:exit, 36.0;
                            description="distance added to cells without phloem under avoidance (cells)"),
                  Parameter(:bias, 0.33;
                            description="change in distance to the richest and poorest quarter of cells under :resource_bias")]
    return CAModel(rule; name=:emerald_ash_borer, states, parameters,
                   description="Mercader et al. (2011) emerald ash borer spread by foraging and phloem")
end

"""
    eab_disperse!(u, p, t; attraction=nothing) -> u

The emergence and dispersal of adults of [Mercader2011](@cite), applied to the
state `u` in place after year `t`, for year `t + 1`. The larvae of one year
(`larvae1`) and the two-year larvae of the matching parity emerge, and the
adults from each cell spread over the cells within `window` cells (clipped at
the edge of the lattice) in proportion to `exp(-b cell (D + y))`, normalised
source by source, where `D` is the distance between cell centres in cells
(`self_distance` for the cell itself) and `y` an extra distance to the
receiving cell:

- `:standard`: none;
- `:avoidance`: `exit` for cells without phloem;
- `:resource_bias` (inferred from the source's description): as avoidance,
  and `D` shortened by `bias` for cells that started with phloem at or above
  the upper quartile over the landscape and lengthened by `bias` for those at
  or below the lower quartile, unless the two quartiles are equal (so a
  uniform landscape has no bias).

If `attraction` is given (a matrix the size of the lattice), each receiving
cell's weight is also multiplied by it before normalising. Adults landing in
border cells are lost.
"""
function eab_disperse!(u, p, t::Integer; attraction=nothing)
    m, n = size(u.phloem)
    odd = isodd(t + 2)   # the source's iteration t + 2
    W = p.window
    D = [hypot(i, j) for i in (-W):W, j in (-W):W]
    D[W + 1, W + 1] = p.self_distance
    extra = zeros(m, n)
    scale = ones(m, n)
    p.behaviour !== :standard && (extra[u.phloem .<= 0] .= p.exit)
    if p.behaviour === :resource_bias
        inside = sort(u.phloem0[.!u.border])
        lo = inside[max(1, ceil(Int, 0.25 * length(inside)))]
        hi = inside[max(1, ceil(Int, 0.75 * length(inside)))]
        for I in eachindex(u.phloem)
            (u.border[I] || hi == lo) && continue
            u.phloem0[I] >= hi && (scale[I] = 1 - p.bias)
            u.phloem0[I] <= lo && (scale[I] = 1 + p.bias)
        end
    end
    adults = zeros(m, n)
    weights = zeros(2W + 1, 2W + 1)
    k = p.b * p.cell
    # without resource bias the weight factorises into the kernel and a factor of the receiving cell
    plain = all(==(1), scale)
    K = exp.(-k .* D)
    F = exp.(-k .* extra)
    attraction === nothing || (F .*= attraction)
    for J in CartesianIndices(u.phloem)
        E = u.larvae1[J] + (odd ? u.larvae2_odd[J] : u.larvae2_even[J])
        E > 0 || continue
        i0, j0 = Tuple(J)
        is = max(1, i0 - W):min(m, i0 + W)
        js = max(1, j0 - W):min(n, j0 + W)
        total = 0.0
        @inbounds for j in js, i in is
            a, c = i - i0 + W + 1, j - j0 + W + 1
            w = plain ? K[a, c] * F[i, j] :
                exp(-k * D[a, c] * scale[i, j]) * F[i, j]
            weights[a, c] = w
            total += w
        end
        @inbounds for j in js, i in is
            adults[i, j] += E * weights[i - i0 + W + 1, j - j0 + W + 1] / total
        end
    end
    adults[u.border] .= 0
    u.adults .= adults
    return u
end

"""
    eab_callback() -> DiscreteCallback

A callback that runs [`eab_disperse!`](@ref) after each year.
"""
function eab_callback()
    affect! = function (integ)
        eab_disperse!(integ.u, integ.p, Int(integ.t))
        return nothing
    end
    return DiscreteCallback((u, t, integ) -> true, affect!; save_positions=(false, false))
end

"""
    eab_environment(kind=:homogeneous; n=205, total=200n^2, rng=Random.Xoshiro(0))
        -> Matrix{Float64}

The synthetic environments of [Mercader2011](@cite), `n × n` cells of phloem
(m²) holding `total` in all: `:homogeneous` (the same in every cell),
`:clumped` (5 × 5 blocks of cells separated by 5 empty cells, one at the
centre), and `:random` (levels 0, 25, ..., 400 m² on equal numbers of cells at
random, then rescaled to the total, with phloem at the centre).
"""
function eab_environment(kind::Symbol=:homogeneous; n::Integer=205,
                         total::Real=200 * n^2, rng=Random.Xoshiro(0))
    isodd(n) || throw(ArgumentError("n must be odd, so the landscape has a centre cell"))
    c = (n + 1) ÷ 2
    A = if kind === :homogeneous
        ones(n, n)
    elseif kind === :clumped
        [mod(i - c + 2, 10) < 5 && mod(j - c + 2, 10) < 5 ? 1.0 : 0.0
         for i in 1:n, j in 1:n]
    elseif kind === :random
        levels = repeat(0:25:400, cld(n^2, 17))[1:(n ^ 2)]
        R = reshape(Float64.(levels[Random.randperm(rng, n^2)]), n, n)
        R[c, c] == 0 && (R[c, c] = 200.0)
        R
    else
        throw(ArgumentError("kind must be :homogeneous, :clumped or :random, got :$kind"))
    end
    return A .* (total / sum(A))
end

"""
    problem_emerald_ash_borer(; environment=:homogeneous, n=205, founders=100,
                              behaviour=:avoidance, years=15, seed=0, pmap=(;),
                              kwargs...) -> CAProblem

The setup of [Mercader2011](@cite): an `environment` of phloem (a name for
[`eab_environment`](@ref), drawn from `seed`, or a matrix) of `n × n` cells of
50 m, with `founders` adults at its centre that emerge and disperse in the
first year, for `years` years. As in the source's code, the landscape is
surrounded by a border of 2 cells of 200 m² of phloem where dispersing
beetles are lost. The model is otherwise deterministic. The state at time `t`
is the landscape after the larvae of year `t` have fed, with the adults that
will lay in year `t + 1` already dispersed.
"""
function problem_emerald_ash_borer(; environment=:homogeneous, n::Integer=205,
                                   founders::Real=100, behaviour::Symbol=:avoidance,
                                   years::Integer=15, seed::Integer=0, pmap=(;), kwargs...)
    behaviour in (:standard, :avoidance, :resource_bias) ||
        throw(ArgumentError("behaviour must be :standard, :avoidance or :resource_bias, got :$behaviour"))
    founders >= 0 || throw(ArgumentError("founders must not be negative, got $founders"))
    A = environment isa Symbol ? eab_environment(environment; n, rng=Random.Xoshiro(seed)) :
        Float64.(environment)
    size(A, 1) == size(A, 2) && isodd(size(A, 1)) ||
        throw(ArgumentError("the environment must be square with a centre cell"))
    n = size(A, 1)
    N = n + 4
    phloem = fill(200.0, N, N)
    phloem[3:(n + 2), 3:(n + 2)] .= A
    border = fill(true, N, N)
    border[3:(n + 2), 3:(n + 2)] .= false
    larvae1 = zeros(N, N)
    larvae1[(N + 1) ÷ 2, (N + 1) ÷ 2] = founders
    model = build_emerald_ash_borer()
    given = merge((; behaviour), (; (Symbol(k) => v for (k, v) in pairs(pmap))...))
    p = merge(defaults(model), given)
    u0 = (; phloem, adults=zeros(N, N), larvae1, larvae2_odd=zeros(N, N),
          larvae2_even=zeros(N, N), border, phloem0=copy(phloem))
    eab_disperse!(u0, p, 0)
    return CAProblem(model, u0, years, given; space=Grid(N, N), seed,
                     callback=eab_callback(), kwargs...)
end

"""
    eab_spread(u; cell=50) -> NamedTuple

The source's measures of an infestation: the `beetles` (all larvae in the
landscape), the `radius` of infestation (the largest distance from the
origin of a cell holding at least one larva, in whole cells times `cell`
metres, as the source's code computes it) and the `median` distance (m)
within which half the larvae lie.
"""
function eab_spread(u; cell::Real=50)
    L = (u.larvae1 .+ u.larvae2_odd .+ u.larvae2_even) .* .!u.border
    N = size(L, 1)
    c = (N + 1) ÷ 2
    beetles = sum(L)
    occupied = findall(>=(1), floor.(L))
    radius = isempty(occupied) ? 0.0 :
             cell * maximum(I -> floor(hypot(I[1] - c, I[2] - c)), occupied)
    d = vec([hypot(I[1] - c, I[2] - c) for I in CartesianIndices(L)])
    order = sortperm(d)
    acc = cumsum(vec(L)[order])
    k = findfirst(>=(beetles / 2), acc)
    median = k === nothing ? 0.0 : cell * d[order][k]
    return (; beetles, radius, median)
end

_register!(:emerald_ash_borer; build=build_emerald_ash_borer,
           problem=problem_emerald_ash_borer)
