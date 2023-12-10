using Test
@time @testset "Acoustics 1" begin
    include("test_acoustics.jl")
end
@time @testset "Acoustics 2" begin
    include("test_acoustics2.jl")
end
@time @testset "Vibroacoustics" begin
    include("test_vibroacoustics.jl")
end
true
