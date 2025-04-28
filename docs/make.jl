using Documenter
using MomentOpt

makedocs(
         sitename="MomentOpt",
         format = Documenter.HTML(
                                  prettyurls = get(ENV, "CI", nothing) == "true"
                                 ),
         pages = [
                  "Introduction" => "index.md",
                  "Getting started" => "started.md",
                  "The Generalized Moment Problem" => "gmp.md",
                  "Certificates of Non Negativity" => "nonneg.md",
                  "A closer look" => "guide.md",
                  "What's in the box?" => "blackbox.md",
                  "API Reference" => "api.md"
                 ],
         modules = [MomentOpt]
        )

deploydocs(
           repo = "github.com/lanl-ansi/MomentOpt.jl.git",
          )
