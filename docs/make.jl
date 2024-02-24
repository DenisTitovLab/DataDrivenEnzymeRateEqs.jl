using EnzymeFitting
using Documenter

DocMeta.setdocmeta!(EnzymeFitting, :DocTestSetup, :(using EnzymeFitting); recursive=true)

makedocs(;
    modules=[EnzymeFitting],
    authors="Denis Titov titov@berkeley.edu and contributors",
    sitename="EnzymeFitting.jl",
    format=Documenter.HTML(;
        canonical="https://DenisTitovLab.github.io/EnzymeFitting.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/DenisTitovLab/EnzymeFitting.jl.git",
    # devbranch="main",
)
