using Documenter
using GyrotronMagnet

makedocs(
    sitename = "GyrotronMagnet",
    authors = "Carlos Hernando",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
    ),
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
    modules = [GyrotronMagnet],
)

deploydocs(
    repo = "github.com/CarlosHernando16/GyrotronMagnet.git",
    devbranch = "main",
)

