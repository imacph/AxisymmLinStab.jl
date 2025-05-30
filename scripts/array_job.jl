using CSV,DataFrames,DelimitedFiles
using AxisymmLinStab
using JSON3
using TimerOutputs

# Create a TimerOutput, this is the main type that keeps track of everything.
const to = TimerOutput()

row_index::Int64 = parse(Int64,ARGS[1])

df = CSV.read("array_config.csv",DataFrame)


parent_directory = @__DIR__
directory = parent_directory*"/folder_"*string(row_index)

n_rad_max::Int64 = df[row_index,:n_rad_max]
n_leg_max::Int64 = df[row_index,:n_leg_max]
rad_ratio::Float64 = df[row_index,:rad_ratio]
ekman::Float64 = df[row_index,:ekman]
freq::Float64 = df[row_index,:freq] 
reynolds::Float64 = df[row_index,:reynolds] 

dt::Float64 = df[row_index,:dt]
n_time_step::Int64 = df[row_index,:n_time_step]

@timeit to "prep_input" (prep_input = AxisymmLinStab.prep_time_step(n_rad_max,n_leg_max,rad_ratio,ekman,freq))


vecs,KE_array = AxisymmLinStab.time_step(dt,n_time_step,reynolds,prep_input,to)

AxisymmLinStab.save_stream_func(vecs*(1+0im),n_rad_max,rad_ratio,n_leg_max,n_leg_max*3,directory)
AxisymmLinStab.save_azim_vel(vecs,n_rad_max,rad_ratio,n_leg_max,n_leg_max*3,directory)

writedlm(directory*"/KE_array.txt",KE_array)



open(directory*"/timeroutput.json", "w") do io
    JSON3.pretty(io, JSON3.write(TimerOutputs.todict(to)))
end