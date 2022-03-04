using DiferentialEvolution
using Documenter

DocMeta.setdocmeta!(DiferentialEvolution, :DocTestSetup, :(using DiferentialEvolution); recursive=true)

makedocs(;
    modules=[DiferentialEvolution],
    authors="Nicolas Urcelay <nico_mx4@hotmail.com> and contributors",
    repo="https://github.com/Urcelay97/DiferentialEvolution.jl/blob/{commit}{path}#{line}",
    sitename="DiferentialEvolution.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Urcelay97.github.io/DiferentialEvolution.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Urcelay97/DiferentialEvolution.jl",
    devbranch="master",
)
