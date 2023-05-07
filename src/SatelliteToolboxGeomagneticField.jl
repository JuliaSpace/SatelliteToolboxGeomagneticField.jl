module SatelliteToolboxGeomagneticField

using LinearAlgebra
using ReferenceFrameRotations
using SatelliteToolboxLegendre
using SatelliteToolboxTransformations
using StaticArrays

############################################################################################
#                                        Constants
############################################################################################

include("./dipole/dipole_coefficients.jl")
include("./igrf/igrf_coefficients.jl")

############################################################################################
#                                         Includes
############################################################################################

include("./dipole/dipole.jl")
include("./igrf/igrf.jl")

end # module SatelliteToolboxGeomagneticField
