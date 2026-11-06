using Asap
using Test

@testset "quadratic fixed-end forces" begin
    # Uniform load. wL/2 and wL²/12 are the fixed-beam end shear and moment.
    L, w = 10.0, 1000.0
    uniform = fem_quadratic_segment(L, 0.0, L, w, 0.0, 0.0)
    @test uniform.R_A ≈ w * L / 2 atol=1e-6
    @test uniform.R_B ≈ w * L / 2 atol=1e-6
    @test uniform.M_A ≈ -w * L^2 / 12 atol=1e-6
    @test uniform.M_B ≈ w * L^2 / 12 atol=1e-6

    # Partial uniform matches the closed-form linear routine.
    partial = fem_quadratic_segment(8.0, 1.0, 5.0, 250.0, 0.0, 0.0)
    closed = Asap._fem_uniform_partial(8.0, 1.0, 5.0, 250.0)
    @test partial.R_A ≈ closed.R_A atol=1e-6
    @test partial.R_B ≈ closed.R_B atol=1e-6
    @test partial.M_A ≈ closed.M_A atol=1e-6
    @test partial.M_B ≈ closed.M_B atol=1e-6

    # Linear ramp w(τ) = 400 τ on [2, 6] of a 10 m member: w0 = 0, w1 = 400.
    ramp = fem_quadratic_segment(10.0, 2.0, 6.0, 0.0, 400.0, 0.0)
    tri = Asap._fem_triangular_partial(10.0, 2.0, 6.0, 400.0)
    @test ramp.R_A ≈ tri.R_A atol=1e-6
    @test ramp.R_B ≈ tri.R_B atol=1e-6
    @test ramp.M_A ≈ tri.M_A atol=1e-6
    @test ramp.M_B ≈ tri.M_B atol=1e-6

    # w(x) = x² on a unit member. Hand integrals:
    # R_A = 1/15, R_B = 4/15, M_A = -1/60, M_B = 1/30.
    quad = fem_quadratic_segment(1.0, 0.0, 1.0, 0.0, 0.0, 1.0)
    @test quad.R_A ≈ 1 / 15 atol=1e-10
    @test quad.R_B ≈ 4 / 15 atol=1e-10
    @test quad.M_A ≈ -1 / 60 atol=1e-10
    @test quad.M_B ≈ 1 / 30 atol=1e-10
    @test quad.R_A + quad.R_B ≈ 1 / 3 atol=1e-10

    # The 12-vector for a zero-curvature piece matches the linear accumulator.
    fem_q = zeros(12)
    fem_l = zeros(12)
    Asap._fem_quadratic_segment!(fem_q, 8.0, 0.25, 0.75, 100.0, 50.0, 0.0, 0.0, 0.0, -1.0)
    Asap._fem_linear_segment!(fem_l, 8.0, 0.25, 0.75, 100.0, 150.0, 0.0, 0.0, -1.0)
    @test fem_q ≈ fem_l atol=1e-6
end
