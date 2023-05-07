# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Tests related to the simplified dipole model.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] http://helios.fmi.fi/~juusolal/geomagnetism/Lectures/Chapter3_dipole.pdf
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# File: ./src/dipole/dipole.jl
# ==========================================================================================

# Function: geomagnetic_dipole_field
# ------------------------------------------------------------------------------------------

@testset "Function geomagnetic_dipole_field" begin
    R0 = 6378.137e3

    # Test 1
    # ======================================================================================
    #
    # Position aligned with the dipole moment vector.

    C = SatelliteToolboxGeomagneticField._GEOMAGNETIC_DIPOLE_MODEL_COEFFICIENTS

    year     = C[32, 1]
    pole_lat = C[32, 2] |> deg2rad
    pole_lon = C[32, 3] |> deg2rad
    m        = C[32, 4] * 1e22

    # Distance from the Earth center.
    r = R0 + 196e3

    # Compute using the simplified equations [1] [nT].
    k0     = 1e-7m
    B_norm = 2k0 / r^3 * 1e9

    # Rotate to ECEF.
    Deg = angle_to_dcm(-(π / 2 - pole_lat), -pole_lon, 0, :YZX)
    B_g = SVector{3}(0, 0, -B_norm)
    B_e_expected = Deg * B_g

    # Compute using the function.
    r_g = SVector{3}(0, 0, r)
    r_e = Deg * r_g

    B_e_result_f64 = geomagnetic_dipole_field(r_e, year)
    B_e_result_f32 = geomagnetic_dipole_field(Float32.(r_e), year)

    @test B_e_expected ≈ B_e_result_f64 atol = 1e-9
    @test eltype(B_e_result_f64) === Float64

    @test B_e_expected ≈ B_e_result_f32 atol = 1e-1
    @test eltype(B_e_result_f32) === Float32

    # Test 2
    # ======================================================================================
    #
    # Position at the magnetic Equator, which must have half the magnitude of that of the
    # test 1.

    B_g = SVector{3}(0, 0, B_norm / 2)
    B_e_expected = Deg * B_g

    r_g = SVector{3}(r, 0, 0)
    r_e = Deg * r_g

    B_e_result_f64 = geomagnetic_dipole_field(r_e, year)
    B_e_result_f32 = geomagnetic_dipole_field(Float32.(r_e), year)

    @test B_e_expected ≈ B_e_result_f64 atol = 1e-9
    @test eltype(B_e_result_f64) === Float64

    @test B_e_expected ≈ B_e_result_f32 atol = 1e-1
    @test eltype(B_e_result_f32) === Float32

    # Test Extrapolation
    # ======================================================================================

    r_e = R0 * (@SVector rand(3))

    B_e_expected = geomagnetic_dipole_field(r_e, 1900)
    B_e_result   = geomagnetic_dipole_field(r_e, 1890)
    @test B_e_result == B_e_expected

    B_e_expected = geomagnetic_dipole_field(r_e, 2026)
    B_e_result   = geomagnetic_dipole_field(r_e, 2025)
    @test B_e_result == B_e_expected
end
