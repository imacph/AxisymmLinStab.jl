using CSV,DataFrames
using AxisymmLinStab

row_index::Int64 = parse(Int64,ARGS[1])

df = CSV.read("array_config.csv",DataFrame)



directory = "folder_"*string(row_index)

n_rad_max::Int64 = df[row_index,:n_rad_max]
n_leg_max::Int64 = df[row_index,:n_leg_max]
rad_ratio::Float64 = df[row_index,:rad_ratio]
ekman::Float64 = df[row_index,:ekman]
freq::Float64 = df[row_index,:freq] 
reynolds::Float64 = df[row_index,:reynolds] 

dt::Float64 = df[row_index,:dt]
n_time_step::Int64 = df[row_index,:n_time_step]

prep_input = AxisymmLinStab.prep_time_step(n_rad_max,n_leg_max,rad_ratio,ekman,freq)


vecs,KE_array = AxisymmLinStab.time_step(dt,n_time_step,reynolds,prep_input)

AxisymmLinStab.save_stream_func(vecs*(1+0im),n_rad_max,rad_ratio,n_leg_max,n_leg_max*3,directory)
AxisymmLinStab.save_azim_vel(vecs,n_rad_max,rad_ratio,n_leg_max,n_leg_max*3,directory)

writedlm(directory*"/KE_array.txt",KE_array)
