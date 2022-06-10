using Distributed
using DataFrames
using CSV
addprocs(10)
begin
    using ProgressMeter
end
df = DataFrame(
        # reporting statistics
            order_p_0=Float64[],
            order_r_0=Float64[],
            order_d_0=Float64[],
            order_p_1=Float64[],
            order_r_1=Float64[],
            order_d_1=Float64[],
            t_relaxation=Float64[],
            propagated_info=Float64[],
        # kinetic parameters
            v=Float64[],
            ℓ=Float64[],
            R=Float64[],
            K=Float64[],
            δ=Float64[],
            s=Float64[],
            n=Int[],
            label=String[],
    )
profile = "uniform_circle_l_v"
overwrite = true
if length(ARGS) == 1
    profile = ARGS[1]
    overwrite = false
elseif length(ARGS) == 2
    profile = ARGS[1]
    overwrite = parse(Bool,ARGS[2])
else
    error("Usage: fish-school-run.jl [profile] [overwrite]")
end
include("profile/$(profile).jl")
# print cmd from cmds
for cmd in cmds
    println(cmd)
end
if isfile("data/$(label).csv") && !overwrite
    println("data found, loading data")
    temp_df = CSV.read("data/$(label).csv",DataFrame)
    for i in 1:length(temp_df.v)
        push!(df,temp_df[i,:])
    end
else
    println("start simulation")
    println(length(cmds))
    @showprogress pmap(run,cmds)
    println("saving data")
    for f in readdir("tmp/")
        temp_df = CSV.read("tmp/$(f)",DataFrame)
        for i in 1:length(temp_df.v)
            push!(df,temp_df[i,:])
        end
    end
    CSV.write("data/$(label).csv",df)
    for f in readdir("tmp/")
        rm("tmp/$(f)")
    end
end
include("fish-school-analysis.jl")
run(`julia fish-school-plot-heatmap.jl $(label) $(field_1) $(field_2)`)
