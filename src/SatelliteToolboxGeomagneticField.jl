module SatelliteToolboxGeomagneticField

using LinearAlgebra
using PrecompileTools
using ReferenceFrameRotations
using SatelliteToolboxLegendre
using SatelliteToolboxTransformations
using StaticArrays

import SatelliteToolboxBase: LowerTriangularStorage, RowMajor

############################################################################################
#                                        Constants                                         #
############################################################################################

include("./dipole/dipole_coefficients.jl")
include("./igrf/igrf_coefficients.jl")

############################################################################################
#                                         Includes                                         #
############################################################################################

include("./dipole/dipole.jl")
include("./igrf/igrf.jl")

############################################################################################
#                                  Precompilation Workload                                  #
############################################################################################

@setup_workload begin
    @compile_workload begin
        P  = Matrix{Float64}(undef, 14, 14)
        dP = similar(P)

        # We use two dates to compile both the interpolation branch (date before the last
        # year with measurements) and the extrapolation branch (date after it).
        for date in (2020.0, 2026.0)
            igrf(date, 6500e3, 0.5, 0.5)
            igrf(date, 400e3, 0.5, 0.5, Val(:geodetic))
            igrf(date, 6500e3, 0.5, 0.5; P = P, dP = dP)
            igrfd(date, 6500e3, 30.0, 30.0)
            igrfd(date, 400e3, 30.0, 30.0, Val(:geodetic))
        end

        geomagnetic_dipole_field(SVector{3, Float64}(7e6, 0.0, 0.0), 2020)
        geomagnetic_dipole_field([7e6, 0.0, 0.0], 2020)
    end
end

end # module SatelliteToolboxGeomagneticField
