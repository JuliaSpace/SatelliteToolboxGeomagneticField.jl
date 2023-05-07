using Test

using DelimitedFiles
using LinearAlgebra
using ReferenceFrameRotations
using SatelliteToolboxGeomagneticField
using StaticArrays

@testset "IGRF" verbose = true begin
    include("./igrf.jl")
end

@testset "Simplified Dipole Model" verbose = true begin
    include("./dipole.jl")
end
