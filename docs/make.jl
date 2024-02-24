using Documenter, FinEtools, FinEtoolsAcoustics

makedocs(
	modules = [FinEtoolsAcoustics],
	doctest = false, clean = true,
	format = Documenter.HTML(prettyurls = false),
	authors = "Petr Krysl",
	sitename = "FinEtoolsAcoustics.jl",
	pages = Any[
	"Home" => "index.md",
	"How to guide" => "guide/guide.md",
	"Types and Functions" => Any[
		"man/man.md"]
		]
	)

deploydocs(
    repo = "github.com/PetrKryslUCSD/FinEtoolsAcoustics.jl.git",
)
