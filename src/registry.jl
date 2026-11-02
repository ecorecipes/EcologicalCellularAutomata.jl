# The model registry. Each model is a source file in `src/models/` that defines a builder (returning
# a `CAModel`) and a problem constructor (the source's default setup), and registers them under its id
# with `_register!`. Its metadata lives in `models/<id>/metadata.toml`, read and checked when the
# package is precompiled; editing a metadata file triggers recompilation through `include_dependency`.

const MODELS_DIR = normpath(joinpath(@__DIR__, "..", "models"))

const CATEGORIES = ("plant disease", "agricultural pest", "forest pest",
                    "bio-geomorphology")
const FIDELITIES = (:reproduction, :adaptation, :inspired)
const REQUIRED_FIELDS = ("id", "title", "category", "organism", "host", "citation", "doi",
                         "bibkey", "ca_role", "fidelity", "validation_in_source",
                         "reproduces", "deviations", "schemes", "source_code")

"""
    ModelInfo

The metadata of a model in the zoo, from its `models/<id>/metadata.toml`:

- `id`, `title`: the model's identifier and a short title.
- `category`: `"plant disease"`, `"agricultural pest"`, `"forest pest"` or
  `"bio-geomorphology"`;
  `organism` and `host`: the pest or pathogen and what it attacks.
- `citation`, `doi`, `bibkey`: the source, checked against the published paper,
  and its key in the package bibliography.
- `ca_role`: what the cellular automaton represents in the source.
- `fidelity`: `:reproduction` (the source's model, with its parameters),
  `:adaptation` (the source's model with stated changes) or `:inspired`
  (built from the source's ideas, not its specification).
- `validation_in_source`: what the source compared its model with.
- `reproduces`: the figure or claim of the source that the model's vignette
  reproduces.
- `deviations`: every way the implementation departs from the source,
  including details the source leaves unstated and how they were filled.
- `schemes`: the update schemes the model is meant to run with.
- `source_code`: the authors' code, if published, and its licence.
- `notes`: anything else.
"""
struct ModelInfo
    id::Symbol
    title::String
    category::String
    organism::String
    host::String
    citation::String
    doi::String
    bibkey::String
    ca_role::String
    fidelity::Symbol
    validation_in_source::String
    reproduces::String
    deviations::Vector{String}
    schemes::Vector{String}
    source_code::String
    notes::String
end

_metadata_error(id, msg) = throw(MetadataError(Symbol(id), msg))

function _read_metadata(path::AbstractString, id::Symbol)
    isfile(path) || _metadata_error(id, "no metadata file at $path")
    d = TOML.parsefile(path)
    for k in REQUIRED_FIELDS
        haskey(d, k) || _metadata_error(id, "the field $k is missing")
    end
    unknown = setdiff(keys(d), (REQUIRED_FIELDS..., "notes"))
    isempty(unknown) || _metadata_error(id, "unknown fields $(sort(collect(unknown)))")
    d["id"] == String(id) || _metadata_error(id, "the file gives the id $(d["id"])")
    d["category"] in CATEGORIES ||
        _metadata_error(id, "the category must be one of $CATEGORIES, got $(d["category"])")
    fidelity = Symbol(d["fidelity"])
    fidelity in FIDELITIES ||
        _metadata_error(id, "the fidelity must be one of $FIDELITIES, got $fidelity")
    for k in ("deviations", "schemes")
        d[k] isa Vector && all(x -> x isa String, d[k]) ||
            _metadata_error(id, "the field $k must be a list of strings")
    end
    startswith(d["doi"], "10.") || _metadata_error(id, "the doi must start with 10.")
    return ModelInfo(id, d["title"], d["category"], d["organism"], d["host"], d["citation"],
                     d["doi"], d["bibkey"], d["ca_role"], fidelity,
                     d["validation_in_source"], d["reproduces"], String.(d["deviations"]),
                     String.(d["schemes"]), d["source_code"], get(d, "notes", ""))
end

struct _Entry
    info::ModelInfo
    build::Function
    problem::Function
end

const REGISTRY = Dict{Symbol,_Entry}()

# Called at the top level of each model file, so the registry is filled during precompilation.
function _register!(id::Symbol; build::Function, problem::Function)
    haskey(REGISTRY, id) && throw(ArgumentError("model $id is registered twice"))
    path = joinpath(MODELS_DIR, String(id), "metadata.toml")
    include_dependency(path)
    REGISTRY[id] = _Entry(_read_metadata(path, id), build, problem)
    return nothing
end

function _entry(id)
    s = Symbol(id)
    haskey(REGISTRY, s) || throw(UnknownModelError(s, list_models()))
    return REGISTRY[s]
end

"""
    list_models(; category=nothing, fidelity=nothing) -> Vector{Symbol}

The ids of the models in the zoo, sorted, optionally only those of one
`category` (`"plant disease"`, `"agricultural pest"`, `"forest pest"`,
`"bio-geomorphology"`) or one
`fidelity` (`:reproduction`, `:adaptation`, `:inspired`).
"""
function list_models(; category=nothing, fidelity=nothing)
    ids = Symbol[]
    for (id, e) in REGISTRY
        category === nothing || e.info.category == category || continue
        fidelity === nothing || e.info.fidelity == Symbol(fidelity) || continue
        push!(ids, id)
    end
    return sort!(ids)
end

"""
    model_info(id) -> ModelInfo

The metadata of the model `id` (a `Symbol` or a string). An unknown id raises
[`UnknownModelError`](@ref).
"""
model_info(id) = _entry(id).info

"""
    build(id; kwargs...) -> CAModel

The `CellularAutomataCore` model `id`, with structural options given as
keywords (documented by each model's builder). The model is backend-agnostic:
solve a problem built from it with any scheme and backend it supports. Its
parameters take the source's values by default and are changed through the
problem's `pmap`.
"""
build(id; kwargs...) = _entry(id).build(; kwargs...)

"""
    problem(id; kwargs...) -> CAProblem

A problem for the model `id` with the setup of its source: the lattice, the
initial state and the run length, which keywords change (documented by each
model's problem constructor). Every model's problem accepts `seed` and `pmap`
(parameter values by name).
"""
problem(id; kwargs...) = _entry(id).problem(; kwargs...)

Base.show(io::IO, m::ModelInfo) = print(io, "ModelInfo(:", m.id, ", \"", m.title, "\")")
function Base.show(io::IO, ::MIME"text/plain", m::ModelInfo)
    println(io, m.id, ": ", m.title)
    println(io, "  ", m.category, ": ", m.organism, " on ", m.host)
    println(io, "  source: ", m.citation, " https://doi.org/", m.doi)
    println(io, "  fidelity: ", m.fidelity)
    println(io, "  reproduces: ", m.reproduces)
    println(io, "  schemes: ", join(m.schemes, ", "))
    print(io, "  deviations: ", length(m.deviations))
    return nothing
end
