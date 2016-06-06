module DOPO
using InformationCapacity
using QHDLJ

export simulate_reservoir, simulate_reservoir2

function simulate_reservoir(data, μ, ν, ρ, Δ, σ; spd=2, hmax=.01, seed=0)
    kappa = 1.
    chi = μ
    beta = 2μ^2/ν
    pump = -ρ/2/ν/sqrt(2)
    println("chi: $chi")
    println("beta: $beta")
    println("Pump amplitude: $pump")

    kappas = [kappa, beta]

    cav = d_opo(kappas, chi, [0, 0.])
    sdata = σ*data
    N = length(data)
    u_t = t->[sdata[round(Int, min(N, max(1, t/Δ)))], pump]
    tlist = linspace(0, N*Δ, N*spd+1)
    z0 = zeros(Complex128, 2)
    nle = solve_nlsystem(cav, z0, tlist, hmax; u_t=u_t, seed=seed)
    us = inputs(nle, ["signal"])
    os = outputs(nle, ["signal", "pump"])
    us, os, nle
end

function simulate_reservoir2(data, χ, β, Δ, σ, u0; spd=2, hmax=.01, seed=0)
    kappa = 1.
    chi = χ
    beta = β
    gamma = (χ^2/β)
    epsmax = 2σ/sqrt(β)*χ
    # max_sig_amp = sqrt((epsmax-kappa/2.)/gamma)
    # max_pump_amp = 2σ/sqrt(β)


    println("γ: $gamma")
    println("epsmax: $epsmax")
    # println("|α_s|_max: $max_sig_amp")
    # println("|α_s|^2_max: $(max_sig_amp^2)")
    # println("|α_p|_max: $max_pump_amp")
    # println("|α_p|^2_max: $(max_pump_amp^2)")



    kappas = [kappa, beta]

    cav = d_opo(kappas, chi, [0, 0.])
    sdata = σ*data
    N = length(data)
    u_t = t->[u0, sdata[round(Int, min(N, max(1, t/Δ)))]]
    tlist = linspace(0, N*Δ, N*spd+1)
    z0 = zeros(Complex128, 2)
    nle = solve_nlsystem(cav, z0, tlist, hmax; u_t=u_t, seed=seed)
    us = inputs(nle, ["pump"])
    os = outputs(nle, ["signal", "pump"])
    us, os, nle
end
end
