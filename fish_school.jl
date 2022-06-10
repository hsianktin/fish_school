using PGFPlotsX
using LinearAlgebra
using CSV
using DataFrames
using Statistics
∑(x) = sum(x)
Λ(x,y) = sum(x'*y) # inner product

function perturb!(ω, S, θ₀)
    for i ∈ S
        ω[2][i] = [cos(θ₀),sin(θ₀)]
    end
end

function state_initiate(n, orientation = "random")
    r̂ = [R*[cos(2π*i/n),sin(2π*i/n)] for i in 1:n]
    if orientation == "random"
        v̂ = [[cos(rand()*2π),sin(rand()*2π)] for i in 1:n]
    elseif orientation == "straight"
        v̂ = [[cos(0),sin(0)] for i in 1:n]
    else
        error("orientation must be random or straight")
    end
    return [r̂, v̂]
end

function normalize(ω)
    return [ω[1],ω[2]/norm(ω[2])]
end

function update(ω, t, ε , dynamics)
    # ω = [r̂, v̂]
    # update the orientation of the fish
    Aₜ = A(ω)
    if dynamics == "vicsek"
        # Aₜ = A(ω) # adjacency matrix
        r̂ = ω[1]
        v̂ = ω[2]
        δr̂ = [[0,0] for i in 1:length(ω[1])] # currently, consider the stationary case
        δv̂ = [[0,0] for i in 1:length(ω[1])]
        # dvᵢ/dt = ∑ⱼ Aᵢⱼ (vⱼ - ⟨vⱼ,vᵢ⟩vᵢ)
        for i = 1:n
            δv̂[i] = ∑(
                [
                    Aₜ[i,j] * (v̂[j] - Λ(v̂[j],v̂[i])*v̂[i] ) for j in 1:n
                ]
            )
        end
        δω = [δr̂, δv̂]
        ωₜ₊₁ = normalize(ω + δω * ε) 
        # print("$(∑(δω))\n")
        return [ωₜ₊₁, t + ε]
    elseif dynamics == "divergence"
        r̂ = ω[1]
        v̂ = ω[2]
        δr̂ = [[0,0] for i in 1:length(ω[1])] # currently, consider the stationary case
        δv̂ = [[0,0] for i in 1:length(ω[1])]
        # dvᵢ/dt = -div(r)(rᵢ - ⟨rᵢ,vᵢ⟩vᵢ) + ∑ⱼAᵢⱼ(vⱼ - ⟨vⱼ,vᵢ⟩vᵢ)
        divₜ = div(ω)
        for i = 1:n
            δv̂[i] = ∑(
                [
                    Aₜ[j,i] * (v̂[j] - Λ(v̂[j],v̂[i])*v̂[i] ) for j in 1:n
                ]
            ) - divₜ * (r̂[i] - Λ(r̂[i], v̂[i])*v̂[i]) * K
        end
        δω = [δr̂, δv̂]
        ωₜ₊₁ = normalize(ω + δω * ε) 
        # print("$(∑(δω))\n")
        return [ωₜ₊₁, t + ε]
    else
        error("dynamics $dynamics not implemented")
    end
end

function A(ω)
    r̂ = ω[1]
    A = zeros(length(r̂), length(r̂))
    # v̂ = ω[2]
    # tg_v̂ = [[cos(ω[i,2]+π/2),sin(ω[i,2]+π/2)] for i in 1:n]
    for i = 1:length(r̂), j = 1:length(r̂)
        if norm(r̂[i]-r̂[j]) < ℓ && i != j
            A[i,j]  = 1
        end
    end
    return A
end
# update law
function div(ω)
    r̂ = ω[1]
    v̂ = ω[2]
    return ∑(Λ.(v̂, r̂))/n
end

function E_p(ω)
    r̂ = ω[1]
    v̂ = ω[2]
    return norm(∑(v̂)/n)
end

function E_r(ω)
    r̂ = ω[1]
    v̂ = ω[2]
    tg_r̂ = [[cos(ω[i,1]+π/2),sin(ω[i,1]+π/2)] for i in 1:n]
    return abs(∑([Λ(tg_r̂[i],v̂[i]) for i in 1:n])/n)
end

n = 20
interval = 100
ε = 1e-3
μ = 10.0
ℓ = 0.4 # absolute length of interaction
R = 1 # radius of group
N = 100
t_f = 10 # only works if plot_flag  == false
K = 1 # measures the effects of divergence
dynamics = "divergence"

if length(ARGS) == 1
    dynamics = ARGS[1]
else
    error("usage: julia fish_school.jl dynamics")
end


for i in 1:N
    global t = 0
    global ω = state_initiate(n, "random")
    # ω₀ = ω
    order_p = Float64[]
    order_r = Float64[]
    order_d = Float64[]
    T = Float64[]
    # adjacency matrix
    global Aₜ = A(ω)
    count = 0
    while t < t_f
        # print("count = $count \r")
        # global count
        # print("$(E_p(ω))\r")
        push!(order_p, E_p(ω))
        push!(order_r, E_r(ω))
        push!(order_d, div_old(ω))
        push!(T, t)
        count += 1
        update!(ω, A, t, ε, μ, dynamics)
    end
    idtfier = rand()
    # print("order_p[end] = $(order_p[end])")
    df = DataFrame(T=[T[i] for i in 1:length(T) if i % interval == 1], order_p=[order_p[i] for i in 1:length(T) if i % interval == 1], order_r=[order_r[i] for i in 1:length(T) if i % interval == 1], order_d=[order_d[i] for i in 1:length(T) if i % interval == 1])
    CSV.write("tmp/sample_$(dynamics)_$(idtfier).csv",df)
    order_p′ = Float64[]
    order_r′ = Float64[]
    order_d′ = Float64[]
    T′ = Float64[]
    θ₀ = 2π * rand()
    S = [i for i in 1:n if abs(ω[i,1] - θ₀) < 2π/n]
    while t < 2t_f 
        push!(order_p′, E_p(ω))
        push!(order_r′, E_r(ω))
        push!(order_d′, div_discrete(ω))
        push!(T′, t-T[end])
        update!(ω, A, t, ε, μ, dynamics)
        perturb!(ω, S, θ₀)
    end
    df′ = DataFrame(T=[T′[i] for i in 1:length(T′) if i % interval == 1], order_p=[order_p′[i] for i in 1:length(T′) if i % interval == 1], order_r=[order_r′[i] for i in 1:length(T′) if i % interval == 1], order_d=[order_d′[i] for i in 1:length(T′) if i % interval == 1])
    CSV.write("tmp/sample_$(dynamics)_$(idtfier)′.csv",df′)
    # print("order_p′[1] = $(order_p′[1])")
    t = 0
    global ω = state_initiate(n, "random")
    order_p′′ = Float64[]
    order_r′′ = Float64[]
    order_d′′ = Float64[]
    T′′ = Float64[]
    θ₀ = 2π * rand()
    S = [i for i in 1:n if abs(ω[i,1] - θ₀) < 2π/n]
    while t < t_f
        push!(order_p′′, E_p(ω))
        push!(order_r′′, E_r(ω))
        push!(order_d′′, div_discrete(ω))
        push!(T′′, t)
        update!(ω, A, t, ε, μ, dynamics)
        perturb!(ω, S, θ₀)
    end
    df′′ = DataFrame(T=[T′′[i] for i in 1:length(T′′) if i % interval == 1], order_p=[order_p′′[i] for i in 1:length(T′′) if i % interval == 1], order_r=[order_r′′[i] for i in 1:length(T′′) if i % interval == 1], order_d=[order_d′′[i] for i in 1:length(T′′) if i % interval == 1])
    CSV.write("tmp/sample_$(dynamics)_$(idtfier)′′.csv",df′′)
end
