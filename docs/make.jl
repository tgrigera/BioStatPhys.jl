push!(LOAD_PATH,"../src/")

using Documenter,BioStatPhys

makedocs(sitename="BioStatPhys.jl",
         authors = "Tomás S. Grigera",
         format = Documenter.HTML( prettyurls = get(ENV, "CI", nothing) == "true" ) ,
         pages=[
             "Home" => "index.md",
             "General statistical tools" => "stat.md"
             ]
         )

if get(ENV, "CI", nothing) == "true" 
    deploydocs( repo = "github.com/tgrigera/BioStatPhys.jl.git")
end

