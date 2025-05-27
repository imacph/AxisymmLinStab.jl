module AxisymmLinStab

using IterativeSolvers
using LinearAlgebra

include("ZerothOrder.jl")

precompile(calc_rad_pts, (Int64, Float64))
precompile(loop_df1_array!, (Array{Float64,2}, Float64, Float64, Int64))
precompile(loop_df2_array!, (Array{Float64,2}, Int64, Vector{Float64}))
precompile(df1_mat, (Int64, Vector{Float64}))
precompile(df2_mat, (Int64, Vector{Float64}))
precompile(calc_diff_mat_array, (Int64, Int64, Vector{Float64}, AbstractArray{Float64}, AbstractArray{Float64}, AbstractArray{Float64}, AbstractArray{Float64}))
precompile(calc_cor_mat_array, (Int64, Int64, Vector{Float64}, AbstractArray{Float64}))
precompile(calc_time_mat_array, (Int64, Int64, Vector{Float64}, AbstractArray{Float64}, AbstractArray{Float64}))
precompile(calc_bc_mat_array, (Int64, Int64))
precompile(mat_array_to_spmat, (AbstractArray{Float64},))
precompile(spatial_mats, (Int64, Int64, Float64))
precompile(order_0_parametrized_system, (Int64, Int64, Float64, Float64, Float64))
precompile(trapezoidal_int, (Vector{Float64}, Int64, Vector{Float64}))

include("SolutionFields.jl")

precompile(calc_theta_grid, (Int64,))
precompile(calc_cart_grid, (Int64, Vector{Float64}, Int64, Vector{Float64}))
precompile(calc_azim_vel, (Vector{ComplexF64}, Int64, Int64, Int64))
precompile(calc_azim_vel, (Vector{Float64}, Int64, Int64, Int64))
precompile(calc_stream_func, (Vector{ComplexF64}, Int64, Int64, Int64))
precompile(calc_stream_func, (Vector{Float64}, Int64, Int64, Int64))
precompile(calc_lsq_stream_func, (Vector{ComplexF64}, Int64, Vector{Float64}, Int64, Int64))
precompile(calc_lsq_stream_func, (Vector{Float64}, Int64, Vector{Float64}, Int64, Int64))
precompile(calc_kinetic_energy, (Vector{Float64}, Int64, Int64, Int64, Float64))
precompile(calc_colat_vel, (Vector{Float64}, Int64, Vector{Float64}, Int64, Int64))
precompile(calc_rad_vel, (Vector{Float64}, Int64, Vector{Float64}, Int64, Int64))
precompile(save_azim_vel, (Vector{ComplexF64}, Int64, Float64, Int64, Int64))
precompile(save_azim_vel, (Vector{Float64}, Int64, Float64, Int64, Int64))
precompile(save_lsq_stream_func, (Vector{ComplexF64}, Int64, Float64, Int64, Int64))
precompile(save_stream_func, (Vector{ComplexF64}, Int64, Float64, Int64, Int64))

include("Coupling.jl")

precompile(Vals_ab_int, (Vector{Float64}, Vector{Float64}, Vector{Float64}, Int64, Vector{Float64}))
precompile(calc_coupling_mat_array, (Int64, Float64, Int64, Int64, Vector{Float64}))
precompile(coupling_matrix, (Int64, Float64, Int64, Int64, Vector{ComplexF64}))
precompile(save_time_average_azim_vel, (Int64, Float64, Int64, Int64, Float64, Float64))


include("time_stepping.jl")

precompile(prep_time_step,(Int64,Int64,Float64,Float64,Float64))
precompile(time_step,(Float64,Int64,Float64,Tuple))   

end
