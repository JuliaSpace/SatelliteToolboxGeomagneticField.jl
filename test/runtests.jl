## Description #############################################################################
#
# Entry point of the test suite.
#
############################################################################################

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

# The quality, performance, and differentiation checks depend on packages that are not
# guaranteed to work on prerelease versions of Julia (JET, AllocCheck, and Zygote).
if isempty(VERSION.prerelease)
    using Aqua
    using AllocCheck
    using JET

    @testset "Performance Tests" verbose = true begin
        include("./performance.jl")
    end

    using ForwardDiff
    using Zygote

    @testset "Zygote Extension" verbose = true begin
        include("./zygote_extension.jl")
    end
else
    @warn "Performance checks and differentiation extension not guaranteed to work on julia-nightly, skipping"
end
