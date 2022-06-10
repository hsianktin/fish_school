
label = "uniform_circle_divergence_delta"
#######################
####### Profile #######
#######################
n = 20 # number of fish
ℓ = 0.4 # absolute length of interaction
R = 1 # radius of group
v = 0.0 # velocity
N = 100
Ks = [exp(i) for i in -1:5] # measures the effects of divergence
δs = [i/10 for i in 0:2:20] # measures noise
s = 0.1 # fraction of stubborn fish in perturbation
IC = "uniform_circle" # initial configuration
dt = 1e-1 # time step
t_f = 2500 # only works if plot_flag  == false
tol = 1e-1 # tolerance for relaxation time
field_1 = "K"
field_2 = "δ"

cmds = []
for K in Ks, δ in δs
    cmd = `julia fish-school-base.jl $n $dt $ℓ $R $v $N $t_f $K $δ $tol $s $IC $label`
    push!(cmds,cmd)
end
