# The registry: every registered model has valid metadata, a builder returning a CAModel and a problem
# constructor; the metadata reader rejects malformed files with MetadataError.

const ECA = EcologicalCellularAutomata

function metadata_text(; drop=nothing, extra="", fields...)
    d = Dict{String,Any}("id" => "toy", "title" => "Toy", "category" => "plant disease",
                         "organism" => "a fungus", "host" => "a plant",
                         "citation" => "A. Author (2000)", "doi" => "10.1000/toy",
                         "bibkey" => "Toy2000", "ca_role" => "contact process",
                         "fidelity" => "reproduction", "validation_in_source" => "none",
                         "reproduces" => "Figure 1", "deviations" => String[],
                         "schemes" => ["Synchronous"], "source_code" => "none")
    for (k, v) in fields
        d[String(k)] = v
    end
    drop === nothing || delete!(d, drop)
    io = IOBuffer()
    ECA.TOML.print(io, d)
    return String(take!(io)) * extra
end

function read_toy(text)
    path = tempname() * ".toml"
    write(path, text)
    return ECA._read_metadata(path, :toy)
end

@testset "registry" begin
    info = read_toy(metadata_text(; deviations=["a", "b"], notes="n"))
    @test info.fidelity === :reproduction && info.deviations == ["a", "b"] &&
          info.notes == "n"
    @test sprint(show, info) == "ModelInfo(:toy, \"Toy\")"
    @test occursin("fidelity: reproduction", sprint(show, MIME"text/plain"(), info))

    @test_throws MetadataError read_toy(metadata_text(; drop="doi"))
    @test_throws MetadataError read_toy(metadata_text(; fidelity="faithful"))
    @test_throws MetadataError read_toy(metadata_text(; category="weed"))
    @test_throws MetadataError read_toy(metadata_text(; id="other"))
    @test_throws MetadataError read_toy(metadata_text(; doi="https://doi.org/10.1/x"))
    @test_throws MetadataError read_toy(metadata_text(; schemes="Synchronous"))
    @test_throws MetadataError read_toy(metadata_text(; extra="colour = \"red\"\n"))
    @test_throws MetadataError ECA._read_metadata(tempname(), :toy)
    @test occursin("model :toy: the field doi is missing",
                   sprint(showerror, MetadataError(:toy, "the field doi is missing")))

    @test_throws UnknownModelError model_info(:no_such_model)
    @test_throws UnknownModelError build("no_such_model")
    @test occursin("no model :x", sprint(showerror, UnknownModelError(:x, [:a, :b])))

    @test issorted(list_models())
    for id in list_models()
        m = model_info(id)
        @test m.id === id && model_info(String(id)) === m
        @test id in list_models(; category=m.category)
        @test id in list_models(; fidelity=m.fidelity)
        @test build(id) isa CAModel
        @test problem(id; seed=1) isa CAProblem
    end
end
