using DataDrivenEnzymeRateEquations
using Documenter

DocMeta.setdocmeta!(DataDrivenEnzymeRateEquations, :DocTestSetup, :(using DataDrivenEnzymeRateEquations); recursive=true)

makedocs(;
    modules=[DataDrivenEnzymeRateEquations],
    authors="Denis Titov titov@berkeley.edu and contributors",
    sitename="DataDrivenEnzymeRateEquations.jl",
    format=Documenter.HTML(;
        canonical="https://DenisTitovLab.github.io/DataDrivenEnzymeRateEquations.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/DenisTitovLab/DataDrivenEnzymeRateEquations.jl.git",
    # devbranch="main",
)
