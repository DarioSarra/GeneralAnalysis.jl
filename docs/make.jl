using GeneralAnalysis
using Documenter

DocMeta.setdocmeta!(GeneralAnalysis, :DocTestSetup, :(using GeneralAnalysis); recursive=true)

makedocs(;
    modules=[GeneralAnalysis],
    authors="Dario",
    repo="https://github.com/DarioSarra/GeneralAnalysis.jl/blob/{commit}{path}#{line}",
    sitename="GeneralAnalysis.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://DarioSarra.github.io/GeneralAnalysis.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/DarioSarra/GeneralAnalysis.jl",
    devbranch="main",
)
