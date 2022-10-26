push!(LOAD_PATH,"../src/")

using Documenter,BioStatPhys

makedocs(sitename="BioStatPhys.jl",
         authors = "Tomás S. Grigera",
         format = Documenter.HTML( prettyurls = get(ENV, "CI", nothing) == "true" ) ,
         pages=[
             "Home" => "index.md",
             "Statistical tools" => "stat.md",
             "Correlation functions" => "corr.md",
             "Miscellaneous tools" => "tools.md"
             ]
         )

if get(ENV, "CI", nothing) == "true" 
    deploydocs( repo = "github.com/tgrigera/BioStatPhys.jl.git")
end
