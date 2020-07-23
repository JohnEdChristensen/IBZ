push!(LOAD_PATH,"../src/")
using Documenter, IBZ

makedocs(sitename="IBZ Documentation",
         format = Documenter.HTML(prettyurls = false),
         modules = [IBZ],
         authors = "Jeremy Jorgensen",
         pages=["index.md", "Documentation.md", "Usage.md"])

# deploydocs(
#    repo = "github.com/jerjorg/IBZ.jl.git",
# )
