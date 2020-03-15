using Documenter
using MomentOpt

makedocs(
         sitename="MomentOpt",
         format = Documenter.HTML(
                                  prettyurls = get(ENV, "CI", nothing) == "true"
                                 ),
         pages = [
                  "Index" => "index.md",
                  "The Generalized Moment Problem" => "gmp.md",
                 ],
         modules = [MomentOpt]
        )

deploydocs(
           repo   = "github.com/lanl-ansi/MomentOpt.jl.git"
          )
