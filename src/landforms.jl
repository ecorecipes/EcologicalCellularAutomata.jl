# Landform geometry: cross-shore profiles built by joining segments end to end. A segment runs
# from one elevation to another over whole cells of 1 m, from the sea to the land; two segments
# join when the first ends at the elevation where the second starts. Joining is associative, and
# the empty segment at an elevation joins without effect, so profiles behave like paths between
# elevations: they compose, and a scenario is a composite of named pieces.

"""
    ProfileSegment(kind, z0, z1, heights)

A piece of a cross-shore profile: the elevations (m) of consecutive 1 m cells,
from the sea to the land, `heights`, starting from elevation `z0` (the seaward
edge) and ending at `z1` (the landward edge, the last cell's elevation for the
built-in pieces). `kind` names the piece (`:shoreface`, `:beach`, `:berm`,
`:dune`, `:hinterland` or any other `Symbol`), which [`reference_profile`](@ref)
reads. Pieces are joined with `*`, in order from the sea, when the first ends
where the second starts; a profile is itself a `ProfileSegment` (of kind
`:profile`) that remembers its pieces.

Built-in pieces: [`shoreface_segment`](@ref), [`slope_segment`](@ref),
[`flat_segment`](@ref) and [`ridge_segment`](@ref).
"""
struct ProfileSegment
    kind::Symbol
    z0::Float64
    z1::Float64
    heights::Vector{Float64}
    pieces::Vector{ProfileSegment}
end
function ProfileSegment(kind::Symbol, z0::Real, z1::Real, heights::AbstractVector{<:Real})
    return ProfileSegment(kind, Float64(z0), Float64(z1), collect(Float64, heights),
                          ProfileSegment[])
end

Base.length(s::ProfileSegment) = length(s.heights)
_pieces(s::ProfileSegment) = isempty(s.pieces) ? [s] : s.pieces
function Base.show(io::IO, s::ProfileSegment)
    kinds = join(map(p -> p.kind, _pieces(s)), " * ")
    return print(io, "ProfileSegment(", kinds, ": ", s.z0, " → ", s.z1, " m over ",
                 length(s),
                 " m)")
end

"""
    a * b -> ProfileSegment

Join two profile pieces, `a` on the seaward side of `b`. The elevation where
`a` ends must be the elevation where `b` starts (to `1e-9` m); otherwise an
`ArgumentError` is thrown. Joining is associative, and an empty piece at an
elevation joins with no effect.
"""
function Base.:*(a::ProfileSegment, b::ProfileSegment)
    isapprox(a.z1, b.z0; atol=1e-9) ||
        throw(ArgumentError("cannot join a piece ending at $(a.z1) m to one starting at $(b.z0) m"))
    isempty(a.heights) && isempty(a.pieces) && return b
    isempty(b.heights) && isempty(b.pieces) && return a
    return ProfileSegment(:profile, a.z0, b.z1, vcat(a.heights, b.heights),
                          vcat(_pieces(a), _pieces(b)))
end

"""
    shoreface_segment(length; A=0.12) -> ProfileSegment

The underwater profile of a beach over `length` m, rising to the shoreline (0
m): a Dean equilibrium profile, depth `A d^(2/3)` at `d` m from the shoreline
(read at the middle of each cell).
"""
function shoreface_segment(length::Integer; A::Real=0.12)
    length >= 0 || throw(ArgumentError("length must be nonnegative, got $length"))
    heights = [-A * (length - k + 0.5)^(2 / 3) for k in 1:length]
    return ProfileSegment(:shoreface, -A * length^(2 / 3), 0.0, heights)
end

"""
    slope_segment(length, z0, z1; kind=:beach) -> ProfileSegment

A straight slope from `z0` to `z1` m over `length` m (the last cell at `z1`).
"""
function slope_segment(length::Integer, z0::Real, z1::Real; kind::Symbol=:beach)
    length >= 0 || throw(ArgumentError("length must be nonnegative, got $length"))
    return ProfileSegment(kind, z0, z1, [z0 + (z1 - z0) * k / length for k in 1:length])
end

"""
    flat_segment(length, z; kind=:berm) -> ProfileSegment

A level stretch at `z` m over `length` m: a berm, a plateau behind a dune
(`kind=:hinterland`), or, with `length = 0`, the empty piece at `z`. Level
pieces at one elevation join into a longer one: `flat_segment(a, z) *
flat_segment(b, z)` has the heights of `flat_segment(a + b, z)`.
"""
function flat_segment(length::Integer, z::Real; kind::Symbol=:berm)
    length >= 0 || throw(ArgumentError("length must be nonnegative, got $length"))
    return ProfileSegment(kind, z, z, fill(Float64(z), length))
end

"""
    ridge_segment(length, foot, crest, back) -> ProfileSegment

A dune ridge over `length` m, rising from `foot` to `crest` m at its middle
and falling to `back` m on its land side, each half a squared sine.
"""
function ridge_segment(length::Integer, foot::Real, crest::Real, back::Real)
    length >= 2 || throw(ArgumentError("a ridge needs at least 2 m, got $length"))
    heights = map(1:length) do k
        s = k / length
        return s < 0.5 ? foot + (crest - foot) * sin(π * s)^2 :
               back + (crest - back) * cos(π * (s - 0.5))^2
    end
    return ProfileSegment(:dune, foot, back, heights)
end

"""
    profile_heights(profile, n; h=nothing) -> Vector

The last `n` cells of a profile (its land end anchored at the last cell), in
metres, or in slabs of height `h` (rounded to `Int32`) if `h` is given. A
profile shorter than `n` is an error.
"""
function profile_heights(p::ProfileSegment, n::Integer; h=nothing)
    length(p) >= n ||
        throw(ArgumentError("the profile has $(length(p)) cells, fewer than $n"))
    z = p.heights[(end - n + 1):end]
    return h === nothing ? z : round.(Int32, z ./ h)
end

"""
    reference_profile(profile) -> ProfileSegment

The surface the sea works towards: the profile with its dunes and hinterland
removed, every piece of kind `:dune` or `:hinterland` replaced by a level
piece at the elevation of the foot of the first dune. It works piece by piece
and keeps the joins: the reference of `a * b` is the reference of `a` joined to
that of `b` when both lie seaward of the dune or both landward.
"""
function reference_profile(p::ProfileSegment)
    ps = _pieces(p)
    k = findfirst(s -> s.kind in (:dune, :hinterland), ps)
    k === nothing && return p
    foot = ps[k].z0
    flat(s) = s.kind in (:dune, :hinterland) ?
              ProfileSegment(s.kind, foot, foot, fill(foot, length(s))) : s
    return foldl(*, map(flat, ps))
end
