using Test

using DelimitedFiles
using LinearAlgebra
using ReferenceFrameRotations
using SatelliteToolboxBase: LowerTriangularStorage, RowMajor
using SatelliteToolboxGeomagneticField
using StaticArrays

@testset "IGRF" verbose = true begin
    include("./igrf.jl")
end

@testset "Simplified Dipole Model" verbose = true begin
    include("./dipole.jl")
end

if isempty(VERSION.prerelease)
    using Pkg
    Pkg.add("JET")
    Pkg.add("AllocCheck")
    Pkg.add("Aqua")

    using JET
    using AllocCheck
    using Aqua

    @testset "Performance Tests" verbose = true begin
        include("./performance.jl")
    end

    Pkg.add("ForwardDiff")
    Pkg.add("Zygote")
    
    using ForwardDiff
    using Zygote

    @testset "Zygote Extension" verbose = true begin
        include("./zygote_extension.jl")
    end
else
    @warn "Performance checks and differentiation extension not guaranteed to work on julia-nightly, skipping"
end
