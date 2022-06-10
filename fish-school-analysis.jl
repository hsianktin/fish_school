# autonomous fish school analysis
# dataframe structure
    # df = DataFrame(
    #         # reporting statistics
    #             order_p_0=Float64[],
    #             order_r_0=Float64[],
    #             order_d_0=Float64[],
    #             order_p_1=Float64[],
    #             order_r_1=Float64[],
    #             order_d_1=Float64[],
    #             t_relaxation=Float64[],
    #             propagated_info=Float64[],
    #         # kinetic parameters
    #             v=Float64[],
    #             ℓ=Float64[],
    #             R=Float64[],
    #             K=Float64[],
    #             δ=Float64[],
    #             s=Float64[],
    #             n=Int[],
    #             label=String[],
    #     )
using DataFrames
using DataFramesMeta
using Statistics
vs = unique(df.v)
ℓs = unique(df.ℓ)
Rs = unique(df.R)
Ks = unique(df.K)
δs = unique(df.δ)
ss = unique(df.s)
ns = unique(df.n)
labels = unique(df.label)

# query the dataframe with given parameters
function query(df,v,ℓ,R,K,δ,s,n,label)
    get_df = @chain df begin
        @rsubset :v == v
        @rsubset :ℓ == ℓ
        @rsubset :R == R
        @rsubset :K == K
        @rsubset :δ == δ
        @rsubset :s == s
        @rsubset :n == n
        @rsubset :label == label
    end
    return get_df
end
function order_class(order_p,order_r)
    if order_p > 0.01 && order_r > 0.01
        return "unknown"
    elseif order_p > 0.01 && order_r ≤ 0.01
        return "polarized"
    elseif order_p ≤ 0.01 && order_r > 0.01
        return "milling"
    elseif order_p ≤ 0.01 && order_r ≤ 0.01
        return "diverging"
    end
end
    
df_report = DataFrame(
    mean_info_propagated = Float64[],
    std_info_propagated = Float64[],
    mean_time_relaxation = Float64[],
    std_time_relaxation = Float64[],
    state_of_order = String[],
    v = Float64[],
    ℓ = Float64[],
    R = Float64[],
    K = Float64[],
    δ = Float64[],
    s = Float64[],
    n = Int[],
    label = String[],
)
if !isfile("data/df_report_$(label).csv") || overwrite
    CSV.write("data/df_report_$(label).csv",df_report)
    for v in vs, ℓ in ℓs, R in Rs, K in Ks, δ in δs, s in ss, n in ns, label in labels
        # query the dataframe with given parameters   
        data = query(df,v,ℓ,R,K,δ,s,n,label)
        # sorting the initial states
        data_filtered = DataFrame(
            state_of_order = String[],
            t_relaxation = Float64[],
            propagated_info = Float64[],
        )
        # fill in data according to the order class
        for j in 1:length(data.order_p_0)
            order_p_0 = data.order_p_0[j]
            order_r_0 = data.order_r_0[j]
            order_d_0 = data.order_d_0[j]
            order_p_1 = data.order_p_1[j]
            order_r_1 = data.order_r_1[j]
            order_d_1 = data.order_d_1[j]
            state_of_order = order_class(order_p_0,order_r_0)
            push!(data_filtered,[
                state_of_order,
                data.t_relaxation[j],
                data.propagated_info[j],
            ])
        end
        for state_of_order in ["unknown", "polarized", "milling", "diverging"]
            # calculate statistics
            data_filtered_filtered = @chain data_filtered begin
                                        @rsubset :state_of_order == state_of_order
                                    end
            data_filtered_filtered_conditioned = @chain data_filtered_filtered begin
                @rsubset :propagated_info > 0.9
            end
            mean_time_relaxation = mean(data_filtered_filtered_conditioned.t_relaxation)
            std_time_relaxation = std(data_filtered_filtered_conditioned.t_relaxation)
            mean_info_propagated = mean(data_filtered_filtered.propagated_info)
            std_info_propagated = std(data_filtered_filtered.propagated_info)
            # fill in dataframe
            push!(df_report,[
                mean_info_propagated,
                std_info_propagated,
                mean_time_relaxation,
                std_time_relaxation,
                state_of_order,
                v,
                ℓ,
                R,
                K,
                δ,
                s,
                n,
                label,
            ])
        end
        # push!(df_report,[
        #     mean(data_filtered.propagated_info),
        #     std(data_filtered.propagated_info),
        #     mean(data_filtered.t_relaxation),
        #     std(data_filtered.t_relaxation),
        #     state_of_order,
        #     v,
        #     ℓ,
        #     R,
        #     K,
        #     δ,
        #     s,
        #     n,
        #     label,
        # ])
    end    
    CSV.write("data/df_report_$(label).csv",df_report;append=true)
else
    println("data found, skipping...")
end
