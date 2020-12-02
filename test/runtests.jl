using GOLFF
using Test
using StaticArrays

@testset "Rotation" begin
    I =  @SMatrix [1.0 0.0 0.0
                   0.0 1.0 0.0
                   0.0 0.0 1.0]

    c = GOLFF.Cartesian(1.0, 2.0, 3.0)
    c′  = GOLFF.rotc(c, I)
    @test c == c′

    m = @SMatrix [0.0 0.0 1.0
                  0.0 1.0 0.0
                  1.0 0.0 0.0]

    c′  = GOLFF.rotc(c, I)
    c′′  = GOLFF.rotc(c′, I)
    @test c == c′′
end