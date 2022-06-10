
label = "uniform_circle_l_divergence"
#######################
####### Profile #######
#######################
n = 20 # number of fish
ℓs = [i/10 for i in 4:2:16] # absolute length of interaction
R = 1 # radius of group
v = 0.0 # velocity
N = 100
Ks = [exp(i) for i in -1:5] # measures the effects of divergence
δ = 0.0 # measures noise
s = 0.1 # fraction of stubborn fish in perturbation
IC = "uniform_circle" # initial configuration
dt = 1e-1 # time step
t_f = 2500 # only works if plot_flag  == false
tol = 1e-1 # tolerance for relaxation time
field_1 = "K"
field_2 = "ℓ"
cmds = []
for K in Ks, ℓ in ℓs
    cmd = `julia fish-school-base.jl $n $dt $ℓ $R $v $N $t_f $K $δ $tol $s $IC $label`
    push!(cmds,cmd)
end

