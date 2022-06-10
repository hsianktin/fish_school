#################
### Profiles ####
#################
profiles = []
# get files ending in .jl and extract their filenames to profiles
for file in readdir("profile")
    if endswith(file,".jl")
        push!(profiles,(file[1:end-3]))
    end
end
print("profiles are:\n")
for profile in profiles
    println("    $(profile) ")
end

##################################
### Compile without overwrite ####
##################################
cmds = []
for profile in profiles
    cmd = `julia fish-school-run.jl $(profile)`
    push!(cmds,cmd)
end
map(run,cmds)
