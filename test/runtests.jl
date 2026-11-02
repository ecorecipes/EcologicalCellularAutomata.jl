using Test
using EcologicalCellularAutomata
using CellularAutomataCore
using Random

@testset "EcologicalCellularAutomata" begin
    include("test_registry.jl")
    include("test_models.jl")
    include("test_errors.jl")
    include("test_docstrings.jl")
end
