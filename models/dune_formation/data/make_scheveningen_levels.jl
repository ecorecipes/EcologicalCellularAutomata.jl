# Builds scheveningen_levels.csv from the measured high and low waters at Scheveningen in Rijkswaterstaat's
# WaterWebservices (no registration). Not run by the package; kept to document how the bundled file was made.
# Needs the JSON package, and fetches one year a second.
#
#   julia make_scheveningen_levels.jl scheveningen_levels.csv
#
# The highest level in each block of 14.6 days from 1 January 1970 (one spring-neap cycle, 25 a year), for
# the first 1250 blocks: 50 years, as Teixeira et al. (2023) describe their forcing. Levels are in metres
# against NAP, the Dutch datum, as measured; no correction is made to mean sea level.
using Dates, Downloads, JSON, Printf

const URL = "https://ddapi20-waterwebservices.rijkswaterstaat.nl/ONLINEWAARNEMINGENSERVICES/OphalenWaarnemingen"
const BLOCK_DAYS = 14.6
const BLOCKS = 1250
const START = DateTime(1970, 1, 1)

function fetch_year(y)
    body = Dict("Locatie" => Dict("Code" => "scheveningen"),
                "AquoPlusWaarnemingMetadata" =>
                    Dict("AquoMetadata" => Dict("Compartiment" => Dict("Code" => "OW"),
                                                "Grootheid" => Dict("Code" => "WATHTE"),
                                                "Groepering" => Dict("Code" => "GETETM2"),
                                                "Hoedanigheid" => Dict("Code" => "NAP"),
                                                "ProcesType" => "meting")),
                "Periode" => Dict("Begindatumtijd" => "$y-01-01T00:00:00.000+01:00",
                                  "Einddatumtijd" => "$y-12-31T23:59:59.000+01:00"))
    out = IOBuffer()
    Downloads.request(URL; method="POST", input=IOBuffer(JSON.json(body)), output=out,
                      headers=["Content-Type" => "application/json",
                               "User-Agent" => "EcologicalCellularAutomata.jl"])
    reply = JSON.parse(String(take!(out)))
    reply["Succesvol"] || error("no data for $y")
    return [m for s in reply["WaarnemingenLijst"] for m in s["MetingenLijst"]]
end

function highest_by_block()
    highest = fill(typemin(Int), BLOCKS)
    last_year = year(START + Millisecond(round(Int64, BLOCKS * BLOCK_DAYS * 86_400_000)))
    for y in year(START):last_year
        for m in fetch_year(y)
            cm = m["Meetwaarde"]["Waarde_Numeriek"]
            (cm === nothing || abs(cm) > 1000) && continue
            stamp = m["Tijdstip"]
            endswith(stamp, "+01:00") || error("unexpected time zone in $stamp")
            t = DateTime(stamp[1:23], dateformat"yyyy-mm-ddTHH:MM:SS.sss")
            k = floor(Int, Dates.value(t - START) / 86_400_000 / BLOCK_DAYS) + 1
            1 <= k <= BLOCKS && (highest[k] = max(highest[k], round(Int, cm)))
        end
        sleep(1)
    end
    any(==(typemin(Int)), highest) && error("a block holds no measurements")
    return highest
end

function main(path)
    highest = highest_by_block()
    open(path, "w") do io
        println(io, "cycle_start,level_m")
        for (k, cm) in enumerate(highest)
            day = Date(START + Millisecond(round(Int64, (k - 1) * BLOCK_DAYS * 86_400_000)))
            println(io, day, ",", @sprintf("%.2f", cm / 100))
        end
    end
    return nothing
end

isempty(ARGS) ? println("usage: julia make_scheveningen_levels.jl output.csv") :
main(ARGS[1])
