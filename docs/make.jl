using DataDrivenEnzymeRateEqs
using Documenter

DocMeta.setdocmeta!(DataDrivenEnzymeRateEqs, :DocTestSetup, :(using DataDrivenEnzymeRateEqs); recursive=true)

makedocs(;
    modules=[DataDrivenEnzymeRateEqs],
    authors="Denis Titov titov@berkeley.edu and contributors",
    sitename="DataDrivenEnzymeRateEqs.jl",
    format=Documenter.HTML(;
        canonical="https://DenisTitovLab.github.io/DataDrivenEnzymeRateEqs.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorial" => "Tutorial.md",
        "API" => "API.md",
    ],
    checkdocs=:exports
)

deploydocs(;
    repo="github.com/DenisTitovLab/DataDrivenEnzymeRateEqs.jl.git",
    # devbranch="main",
)
