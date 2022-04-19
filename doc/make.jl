using Documenter, pomin

makedocs(
    sitename    =   "pomin.jl",
    format      = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    authors     = "Justin C. Feng and Mark Baumann",
    pages = [
                "Home"              => "index.md",
                "Hamlsb"            => "hamlsb.md",
                "Integrators"       => "integrators.md",
            ]
        )
