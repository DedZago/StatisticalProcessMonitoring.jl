push!(LOAD_PATH, "../src/")
using StatisticalProcessMonitoring
using Documenter
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "src", "biblio.bib"))


DocMeta.setdocmeta!(StatisticalProcessMonitoring, :DocTestSetup, :(using StatisticalProcessMonitoring); recursive=true)

makedocs(;
   modules=[StatisticalProcessMonitoring],
   authors="Daniele Zago <daniele.zago.1@phd.unipd.it>",
   repo="https://github.com/DedZago/StatisticalProcessMonitoring.jl/blob/{commit}{path}#{line}",
   sitename="StatisticalProcessMonitoring.jl",
   plugins=[bib],
   format=Documenter.HTML(;
                          prettyurls=get(ENV, "CI", "false") == "true",
                          canonical="https://DedZago.github.io/StatisticalProcessMonitoring.jl",
                          edit_link="main",
                          assets=String[],
                         ),
   pages=[
          "Introduction" => "index.md",
          "Theory" => "theory.md",
          "Getting started" => "getting-started.md",
          # "Basic usage" => ["using_control_charts.md"],
          "Examples" => ["monitoring_mean_covariance.md", "monitoring_autoregressive.md", "monitoring_risk_adjusted.md", "monitoring_nonparametric_profiles.md"],
          "Reference" => ["statistics.md", "control_limits.md", "nominal_properties.md", "phase_2.md", "control_charts.md","optimization.md"]
          "Bibliography" => ["bibliography.md"]
         ],
  )

deploydocs(;
     repo="github.com/DedZago/StatisticalProcessMonitoring.jl",
     devbranch="main",
    )
