# This file needs to be put into the local installation of SPICE

include("../../src/SPICE.jl")
using .SPICE

# ---- State-change matrix
const ν = Array([1 -1]')
# ---- Initial state-space
x0 = [0]

# ---- Bounds on the initial search space for parameters
const bounds = [[1e-5, 1e-8], [1000, 10]]

# ---- Highest orders of reaction for tau-leaping
const hor = [1]

# ---- Hazard functions
function F(p)
	H(p)
    p.a .*= p.θ
end

function H(p)
    p.a[1] = 1
    p.a[2] = p.x[1]
end

# ---- Configure system in SPICE
system = System(Model(x0,F,H,ν,bounds, hor=hor, obsname=[:X]), raw"path/to/where/source.py/saves/simulation/CSV", routine=CEM(ssa=:Tau, nElite = 10, nRepeat = 1, nSamples=1000, maxIter=250, mSamples=20000, shoot=false, splitting=false, tauOpt=TauOpt(ϵ=0.1)))

# ---- Call parameter estimation routine
estimate(system, 1, "srate100_degrate1_test")