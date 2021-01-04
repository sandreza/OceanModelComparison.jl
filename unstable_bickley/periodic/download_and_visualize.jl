function remote_download(server, directory, file)
    commandline = server * directory * file
    Base.run(`rsync -av $commandline .`)
    return nothing
end
server = "sandre@tartarus.mit.edu:~"
directory = "/Desktop/Julia/OceanModelComparison.jl/" 
file = "climate_machine_unstable_bickley_jet_Ne43_Np2_ν0.0e+00_no_rotation_diffusive_cfl1.0e-02_reduction1.0e+04_smoothness_exponent2.0e+00.jld2"
file = "climate_machine_unstable_bickley_jet_Ne43_Np2_ν0.0e+00_no_rotation_diffusive_cfl1.0e-02_reduction1.0e+00_smoothness_exponent2.0e+00.jld2"
file = "climate_machine_unstable_bickley_jet_Ne43_Np2_ν0.0e+00_no_rotation_diffusive_cfl1.0e-03_reduction1.0e+01_smoothness_exponent1.0e+00.jld2"
file = "climate_machine_unstable_bickley_jet_Ne32_Np3_ν0.0e+00_no_rotation_diffusive_cfl1.0e-01_reduction1.0e+04_smoothness_exponent1.0e+00.jld2"
file = "climate_machine_unstable_bickley_jet_Ne21_Np5_ν0.0e+00_no_rotation_diffusive_cfl1.0e-01_reduction1.0e+01_smoothness_exponent2.0e+00.jld2"
file = "climate_machine_unstable_bickley_jet_Ne18_Np6_ν0.0e+00_no_rotation_diffusive_cfl1.0e+00_reduction1.0e+04_smoothness_exponent1.0e+00.jld2"
file = "climate_machine_unstable_bickley_jet_Ne18_Np6_ν0.0e+00_no_rotation_diffusive_cfl1.0e+00_reduction1.0e+08_smoothness_exponent1.0e+00.jld2"
file = "climate_machine_unstable_bickley_jet_Ne73_Np6_ν0.0e+00_no_rotation_diffusive_cfl1.0e+00_reduction1.0e+01_smoothness_exponent1.0e+00.jld2"
file = "climate_machine_unstable_bickley_jet_Ne37_Np6_ν0.0e+00_no_rotation_diffusive_cfl1.0e+00_reduction1.0e+08_smoothness_exponent1.0e+00.jld2"
file = "climate_machine_unstable_bickley_jet_Ne37_Np6_ν0.0e+00_no_rotation_diffusive_cfl1.0e+00_reduction1.0e+01_smoothness_exponent1.0e+00.jld2"
file = "climate_machine_unstable_bickley_jet_Ne73_Np6_ν0.0e+00_no_rotation_diffusive_cfl1.0e+00_reduction1.0e+08_smoothness_exponent1.0e+00.jld2"
file = "climate_machine_unstable_bickley_jet_Ne64_Np3_ν0.0e+00_no_rotation_diffusive_cfl1.0e+00_reduction1.0e+08_smoothness_exponent2.0e+00.jld2"
file = "climate_machine_unstable_bickley_jet_Ne128_Np3_ν0.0e+00_no_rotation_diffusive_cfl1.0e+00_reduction1.0e+08_smoothness_exponent2.0e+00.jld2"
file = "climate_machine_unstable_bickley_jet_Ne51_Np4_ν0.0e+00_no_rotation_diffusive_cfl1.0e+00_reduction1.0e+08_smoothness_exponent1.0e+00.jld2"
file = "climate_machine_unstable_bickley_jet_Ne102_Np4_ν0.0e+00_no_rotation_diffusive_cfl1.0e+00_reduction1.0e+08_smoothness_exponent1.0e+00.jld2"
file = "climate_machine_unstable_bickley_jet_Ne51_Np4_ν0.0e+00_no_rotation_diffusive_cfl1.0e-01_reduction1.0e+01_smoothness_exponent2.0e+00.jld2"
remote_download(server, directory, file)
##
DOF = 512 
Np = 6
Ne = round(Int, DOF / (Np+1))
effective_node_spacing(Ne, Np, Lx=4π) = Lx / (Ne * (Np + 1)^2)
time_step = 1.0 * effective_node_spacing(Ne, Np) / c 
κ = effective_node_spacing(Ne, Np)^2 / time_step
##
visualize(file[1:end-5])

