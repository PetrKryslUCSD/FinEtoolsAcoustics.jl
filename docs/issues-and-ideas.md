
Issues and ideas:


-- Documenter:
using FinEtoolsAcoustics
using DocumenterTools
Travis.genkeys(user="PetrKryslUCSD", repo="https://github.com/PetrKryslUCSD/FinEtoolsAcoustics.jl")

using Pkg; Pkg.add("DocumenterTools");                                 
using DocumenterTools                                                  
DocumenterTools.genkeys(user="PetrKryslUCSD", repo="git@github.com:PetrKryslUCSD/FinEtoolsAcoustics.jl.git")                                                  
using Pkg; Pkg.rm("DocumenterTools");  