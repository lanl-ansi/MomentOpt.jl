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
                  "What's in the box?" => "blackbox.md"
                 ],
         modules = [MomentOpt]
        )

deploydocs(
           repo   = "github.com/lanl-ansi/MomentOpt.jl.git"
          )
