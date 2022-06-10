using LinearAlgebra
using CSV
using DataFrames
using Statistics

##########################
#### Utils ###############
##########################

# sum function that preserve element types
function ∑(x)
    # this is a hack to preserve the type of the elements
    ∑x = x[1]
    for i in 2:length(x)
        ∑x = ∑x + x[i]
    end
    return ∑x
end

# normailze the velocity vector of the state
function normalize(ω)
    return [ω[1], [ω[2][i] / norm(ω[2][i]) for i in 1:length(ω[2])]]
end

Λ(x, y) = sum(x' * y) # inner product

function perturb(ω, S, θ₀)
    η = ω
    for i ∈ S
        η[2][i] = [cos(θ₀), sin(θ₀)]
    end
    return η
end

# generate a random velocity vector
function V() # random 
    θ = rand() * 2π
    return [cos(θ), sin(θ)]
end
# get a vector according to the angle θ
function V(θ)
    return [cos(θ), sin(θ)]
end


# ∂ₜω
function ∂ₜω(ω, dynamics, δ, K, ℓ)
    Aₜ = A(ω, ℓ)
    r̂ = ω[1]
    v̂ = ω[2]
    δr̂ = v̂ * v # currently, consider the stationary case
    δv̂ = [[0.0, 0.0] for i in 1:length(ω[1])]
    if dynamics == "vicsek"
        # dvᵢ/dt = ∑ⱼ Aᵢⱼ (vⱼ - ⟨vⱼ,vᵢ⟩vᵢ)
        for i = 1:n
            δv̂[i] = ∑(
                [
                Aₜ[i, j] * (v̂[j] - Λ(v̂[j], v̂[i]) * v̂[i]) for j in 1:n
            ]
            ) + δ * V()
        end
        ∂ₜω = [δr̂, δv̂]
        return ∂ₜω
    elseif dynamics == "divergence"
        # dvᵢ/dt = -div(r)(rᵢ - ⟨rᵢ,vᵢ⟩vᵢ) + ∑ⱼAᵢⱼ(vⱼ - ⟨vⱼ,vᵢ⟩vᵢ) + δv
        divₜ = div(ω)
        for i = 1:n
            δv̂[i] = ∑(
                [
                Aₜ[j, i] * (v̂[j] - Λ(v̂[j], v̂[i]) * v̂[i]) for j in 1:n
            ]
            ) - K * divₜ * (
                        r̂[i] - Λ(r̂[i], v̂[i]) * v̂[i]
                    ) + δ * V()
        end
        ∂ₜω = [δr̂, δv̂]
        return ∂ₜω
    else
        error("dynamics $dynamics not implemented")
    end
end

# empoly the reflecting boundary condition
function mod_2R(ω)
    local n = length(ω[1])
    for i in 1:n
        while ω[1][i][1] > 2R
            ω[1][i][1] -= 2R
        end
        while ω[1][i][1] < 2R
            ω[1][i][1] += 2R
        end
        while ω[1][i][2] > 2R
            ω[1][i][2] -= 2R
        end
        while ω[1][i][2] < 2R
            ω[1][i][2] += 2R
        end
    end
    return ω
end

# one-step integration ∂ₜωdt
function ∂ₜωdt(ω, t, dt, dynamics, type, δ, K, ℓ)
    # ω = [r̂, v̂]
    # update the orientation of the fish
    if type == "normal"
        ωₜ₊₁ = normalize(ω + ∂ₜω(ω, dynamics, δ, K, ℓ) * dt)
        # if v > 0
        #     ωₜ₊₁ = mod_2R(ωₜ₊₁)
        # end
        return [ωₜ₊₁, t + dt]
    elseif type == "perturb"
        ωₜ₊₁ = normalize(ω + ∂ₜω(ω, dynamics, δ, K, ℓ) * dt)
        # if v > 0
        #     ωₜ₊₁ = mod_2R(ωₜ₊₁)
        # end
        return [perturb(ωₜ₊₁, S, θ₀), t + dt]
    end
end

# integration of ∂ₜωdt until converge
function ∫∂ₜωdt(ω, t, dt, t_f, dynamics, type, δ, K, ℓ, tol)
    Ω = [] # record the average orientation
    order_p = Float64[]
    order_r = Float64[]
    order_d = Float64[]
    T = Float64[]
    count = 0
    while t < 1
        # print("count = $count \r")
        # count
        # print("$(E_p(ω))\r")
        push!(Ω, ω)
        push!(order_p, E_p(ω))
        push!(order_r, E_r(ω))
        push!(order_d, div(ω))
        push!(T, t)
        count += 1
        ω, t = ∂ₜωdt(ω, t, dt, dynamics, type, δ, K, ℓ)
    end
    while (t < t_f && norm(∂ₜωdt(ω, t, dt, dynamics, type, δ, K, ℓ)[1][2] .- ω[2]) > 1e-10 * dt) || t < 10
        # print("count = $count \r")
        # count
        # print("$(E_p(ω))\r")
        push!(Ω, ω)
        push!(order_p, E_p(ω))
        push!(order_r, E_r(ω))
        push!(order_d, div(ω))
        push!(T, t)
        count += 1
        ω, t = ∂ₜωdt(ω, t, dt, dynamics, type, δ, K, ℓ)
    end
    # println(t)
    return [
        ω,
        order_p[1],
        order_r[1],
        order_d[1],
        order_p[end],
        order_r[end],
        order_d[end],
        relaxation_time(T, Ω, tol)
    ]
end

# get adjacency matrix
function A(ω, ℓ)
    r̂ = ω[1]
    A = zeros(length(r̂), length(r̂))
    for i = 1:length(r̂), j = 1:length(r̂)
        if norm(r̂[i] - r̂[j]) < ℓ && i != j
            A[i, j] = 1
        end
    end
    return A
end

# get ``divergence''
function div(ω)
    r̂ = ω[1]
    v̂ = ω[2]
    return ∑(Λ.(v̂, r̂)) / n
end

# order parameter of polarity
function E_p(ω)
    r̂ = ω[1]
    v̂ = ω[2]
    return norm(∑(v̂) / n)
end

# order parameter of rotation
function E_r(ω)
    r̂ = ω[1]
    v̂ = ω[2]
    tg_r̂ = [[-r̂[i][2], r̂[i][1]] for i in 1:n]
    return abs(∑([Λ(tg_r̂[i], v̂[i]) for i in 1:n]) / n)
end

# compute the time till relaxation
function relaxation_time(T, Ω, tol=1e-1)
    # based on trajectories of order_p and order_r as a function of t ∈ T, determine the relaxation time
    ω_f = Ω[end]
    diffs = [norm(Ω[i][1] .- ω_f[1]) + norm(Ω[i][2] .- ω_f[2]) for i in 1:length(Ω)]./(norm(ω_f[1]) + norm(ω_f[2]))
    i = 1
    t0 = T[1]
    t = T[end]
    while diffs[end-i] < tol && length(diffs) - i > 1
        i += 1
        t = T[end-i]
    end
    return t - t0
end



function main(
    dynamics::String, # dynamics
    n::Int, # number of fish
    dt::Number, # time step
    ℓ::Number, # absolute length of interaction
    R::Number, # radius of group
    v::Number, # velocity
    N::Number, # number of samples
    t_f::Number, # maximum time for integration
    K::Number, # measure of the effects of divergence
    δ::Number, # measure of noise
    tol::Number, # tolerance for relaxation time
    IC::String, # function for initial condition
    label::String, # label for profile
)
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
    
    for i in 1:N
        t = 0
        type = "normal"
        ω₀ = initialCondition(n, R, IC)
        ω₁, order_p_0, order_r_0, order_d_0, order_p_1, order_r_1, order_d_1, relax_time = ∫∂ₜωdt(ω₀, t, dt, t_f, dynamics, type, δ, K, ℓ, tol)
        # push!(df,[
        #     order_p_0,
        #     order_r_0,
        #     order_d_0,
        #     order_p_1,
        #     order_r_1,
        #     order_d_1,
        #     relax_time,
        #     type,
        #     v,
        #     ℓ,
        #     R,
        #     K,
        #     δ,
        #     s,
        #     n,
        #     label])
        type = "perturb"
        # additional parameter for perturbation is broadcast to global scope
        t = 0 # reset time
        ω₂, order_p_0, order_r_0, order_d_0, order_p_1, order_r_1, order_d_1, relax_time = ∫∂ₜωdt(ω₁, t, dt, t_f, dynamics, type, δ, K, ℓ, tol)
        push!(df, [
            order_p_0,
            order_r_0,
            order_d_0,
            order_p_1,
            order_r_1,
            order_d_1,
            relax_time,
            Λ(∑(ω₂[2]) / n, V(θ₀)),
            v,
            ℓ,
            R,
            K,
            δ,
            s,
            n,
            label])
    end
    return df
end



if length(ARGS) == 13
    n = parse(Int, ARGS[1])
    dt = parse(Float64, ARGS[2])
    ℓ = parse(Float64, ARGS[3])
    R = parse(Float64, ARGS[4])
    v = parse(Float64, ARGS[5])
    N = parse(Int, ARGS[6])
    t_f = parse(Float64, ARGS[7])
    K = parse(Float64, ARGS[8])
    δ = parse(Float64, ARGS[9])
    tol = parse(Float64, ARGS[10])
    s = parse(Float64, ARGS[11])
    IC = ARGS[12]
    label = ARGS[13]
else
    error("usage: julia fish-school-base.jl \$n \$dt \$ℓ \$R \$v \$N \$t_f \$K \$δ \$tol \$s \$IC \$label")
end
# hard coding ``uninteresting'' parameters and functions
dynamics = "divergence" # vicsek is subsumed by divergence with K = 0
S = [i for i in 1:floor(Int, n * s)]
θ₀ = θ = rand() * 2π
# initial condition
function initialCondition(n, R=1, IC="uniform_circle")
    if IC == "uniform_circle"
        r̂ = [R * [cos(2π * i / n), sin(2π * i / n)] for i in 1:n]
        v̂ = [V() for i in 1:n]
    elseif IC == "random_box"
        r̂ = [R * [rand(), rand()] for i in 1:n]
        v̂ = [V() for i in 1:n]
    else
        error("IC $(IC) not implemented")
    end
    return [r̂, v̂]
end

CSV.write("./tmp/$(label)_$(rand()).csv", main(dynamics, n, dt, ℓ, R, v, N, t_f, K, δ, tol, IC, label))
