push!(LOAD_PATH, "../src/")
using SPM
using Documenter

DocMeta.setdocmeta!(SPM, :DocTestSetup, :(using SPM); recursive=true)

makedocs(;
    modules=[SPM],
    authors="Daniele Zago <daniele.zago.1@phd.unipd.it>",
    repo="https://github.com/DedZago/SPM.jl/blob/{commit}{path}#{line}",
    sitename="SPM.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://DedZago.github.io/SPM.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/DedZago/SPM.jl",
    devbranch="main",
)
