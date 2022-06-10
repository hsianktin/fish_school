using PGFPlotsX
using DataFrames
using DataFramesMeta
using Statistics
using CSV
include("fish-school-utils.jl")
fields = ["K","ℓ"]
label = "uniform_circle"
if length(ARGS) == 1
    label = ARGS[1]
elseif length(ARGS) == 3
    label = ARGS[1]
    field_1 = ARGS[2]
    field_2 = ARGS[3]
    fields = [field_1,field_2]
else
    error("profile name must be specified")
end

df = CSV.read("data/df_report_$(label).csv",DataFrame)
for col in eachcol(df)
    replace!(col,NaN => 0)
end
df_plot = @chain df begin
    @rsubset :label == label
end
sort!(df_plot,[Symbol(fields[1]),Symbol(fields[2])])
CSV.write("fig/df_plot_$(label).csv",df_plot)
states = unique(df_plot.state_of_order)
push!(PGFPlotsX.CUSTOM_PREAMBLE,raw"\usepgfplotslibrary{colormaps}")
push!(PGFPlotsX.CUSTOM_PREAMBLE,raw"\usetikzlibrary[pgfplots.colormaps]")
# mean information propagated
for state in states
    df_plot_current_state = @chain df_plot begin
        @rsubset :state_of_order == state
    end
    CSV.write("fig/df_plot_$(label)_$(state).csv",df_plot_current_state)
    if unique(df_plot_current_state.mean_info_propagated) == [0]
        # refuse to plot because no meanningful data
        continue
    end
    t = @pgf Table({"col sep"="comma"},"df_plot_$(label)_$(state).csv")
    t["x"] = fields[1]
    t["y"] = fields[2]
    t["z"] = "mean_info_propagated"
    # t["discard if not={state_of_order}{$(state)}"] = nothing
    # select certain parameter in PGFPlotsX


    p = @pgf SemiLogXAxis({
        width = "3.4in",
        height = "3.4in",
        xlabel = unicode2latex(fields[1]),
        ylabel = unicode2latex(fields[2]),
        title = "mean information propagated",
        view = (0, 90),
        colorbar,
        "point meta min"= 0,
        "point meta max"= 1,
        "colormap/temp",
        "clip"="false",
    },
        Plot3({
            surf,
            shader = "flat",
            "mesh/rows" = length(unique(df_plot_current_state[:,Symbol(fields[1])])),
            "mesh/cols" = length(unique(df_plot_current_state[:,Symbol(fields[2])])),
        },t),
        Legend(state)
    )
    pgfsave("fig/plot_mean_info_propagated_$(label)_$(state).tex",p)
end

# mean relaxation time
for state in states
    df_plot_current_state = @chain df_plot begin
        @rsubset :state_of_order == state
    end
    CSV.write("fig/df_plot_$(label)_$(state).csv",df_plot_current_state)
    if unique(df_plot_current_state.mean_time_relaxation) == [0]
        # refuse to plot because no meanningful data
        continue
    end
    t = @pgf Table({"col sep"="comma"},"df_plot_$(label)_$(state).csv")
    t["x"] = fields[1]
    t["y"] = fields[2]
    t["z"] = "mean_time_relaxation"
    # t["discard if not={state_of_order}{$(state)}"] = nothing
    # select certain parameter in PGFPlotsX


    p = @pgf SemiLogXAxis({
        "restrict y to domain"="0.4:1.6",
        width = "3.4in",
        height = "3.4in",
        xlabel = unicode2latex(fields[1]),
        ylabel = unicode2latex(fields[2]),
        title = "mean time relaxation",
        view = (0, 90),
        colorbar,
        "colormap/jet",
        "clip"="false",
    },
        Plot3({
            surf,
            shader = "flat",
            "mesh/rows" = length(unique(df_plot_current_state[:,Symbol(fields[1])])),
            "mesh/cols" = length(unique(df_plot_current_state[:,Symbol(fields[2])])),
        },t),
        Legend(state)
    )
    pgfsave("fig/plot_mean_time_relaxation_$(label)_$(state).tex",p)
end
