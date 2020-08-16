using Literate

for t in readdir(".")
    if occursin(r".*_tut.jl", t)
        println("\nTutorial $t in $(pwd())\n")
        Literate.markdown(t, "."; documenter=false);
        # Literate.notebook(t, "."; execute=false, documenter=false);
        # @show n, ext = splitext(t)
        # @show "../docs/tutorials/" * n * ".md"
        # cp(n * ".md", "../docs/src/tutorials/" * n * ".md", force = true)
    end
end