using Test

using AllocCheck
using Aqua
using DelimitedFiles
using ForwardDiff
using JET
using LinearAlgebra
using ReferenceFrameRotations
using SatelliteToolboxBase: LowerTriangularStorage, RowMajor
using SatelliteToolboxGeomagneticField
using StaticArrays
using Zygote

@testset "IGRF" verbose = true begin
    include("./igrf.jl")
end

@testset "Simplified Dipole Model" verbose = true begin
    include("./dipole.jl")
end

@testset "Performance" verbose = true begin
    include("./performance.jl")
end

@testset "Zygote Extension" verbose = true begin
    include("./zygote_extension.jl")
end
