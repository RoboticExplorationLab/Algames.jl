using Documenter
using Algames  # your package name here

makedocs(
    sitename = "Algames",  # your package name here
    format = Documenter.HTML(prettyurls = false),  # optional
    pages = [
        "Introduction" => "index.md",
        "Basics" => [
            "Creating a New Package" => "new_repo.md",
            "Adding Dependencies" => "deps.md",
            "Adding Tests" => "tests.md",
            "Setting up CI" => "CI.md",
            "Adding Code Coverage" => "codecov.md",
            "Building documentation" => "docs.md",
            "Registering a Package" => "registering.md"
        ],
        "Sample Documentation" => "package_docs.md"
    ]
)

deploydocs(
    repo = "github.com/simon-lc/Algames.jl.git",
)
