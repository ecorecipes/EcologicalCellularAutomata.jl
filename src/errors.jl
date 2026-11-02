# Exceptions of EcologicalCellularAutomata. Every exception the package defines lives in this
# file, with its `showerror` method, and subtypes the root `EcologicalCellularAutomataError`.
#
# - Carry the offending names (layer, rule, backend, model) in the fields and the message.
#   Library code never calls a bare `error("...")`.
# - Invalid arguments and keywords raise `ArgumentError`, not a type from this file.
# - A docstring names another package's type as a code span, never with `@ref`.

"""
    EcologicalCellularAutomataError

Abstract supertype of every exception that EcologicalCellularAutomata defines. It marks the
package that introduced an error, not a kind of failure. Invalid arguments and
keywords raise Base's `ArgumentError` instead, which is outside this root.
"""
abstract type EcologicalCellularAutomataError <: Exception end

"""
    UnknownModelError(id, known)

Raised when a model id is not in the zoo; `known` lists the ids that are.
"""
struct UnknownModelError <: EcologicalCellularAutomataError
    id::Symbol
    known::Vector{Symbol}
end
function Base.showerror(io::IO, e::UnknownModelError)
    return print(io, "UnknownModelError: the zoo has no model :", e.id, "; models are ",
                 join((":" * String(k) for k in e.known), ", "))
end

"""
    MetadataError(id, message)

Raised, when the package is compiled, if the metadata file of the model `id`
is missing, lacks a required field, or has a field with an invalid value.
"""
struct MetadataError <: EcologicalCellularAutomataError
    id::Symbol
    message::String
end
function Base.showerror(io::IO, e::MetadataError)
    return print(io, "MetadataError: model :", e.id, ": ", e.message)
end
