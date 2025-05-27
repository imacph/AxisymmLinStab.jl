using AxisymmLinStab
using Test

@testset "AxisymmLinStab.jl" begin
    @test 1 + 1 == 2

    # Test calc_rad_pts
    r_i, r_o, rad_pts = AxisymmLinStab.calc_rad_pts(5, 0.5)
    @test length(rad_pts) == 5
    @test r_i < r_o
    @test all(isfinite, rad_pts)

    # Test df1_mat and df2_mat
    df1 = AxisymmLinStab.df1_mat(5, rad_pts)
    df2 = AxisymmLinStab.df2_mat(5, rad_pts)
    @test size(df1, 1) == 5 && size(df1, 2) == 5
    @test size(df2, 1) == 5 && size(df2, 2) == 5

    # Test spatial_mats returns expected number of matrices
    mats = AxisymmLinStab.spatial_mats(5, 5, 0.5)
    @test length(mats) == 4

    # Test order_0_parametrized_system returns expected types
    matsys = AxisymmLinStab.order_0_parametrized_system(5, 5, 0.5, 1e-3, 2.0)
    @test length(matsys) == 2 

    # Test calc_theta_grid
    theta = AxisymmLinStab.calc_theta_grid(6)
    @test length(theta) == 6
    @test all(x -> x > 0 && x < π, theta)


end
