## Description #############################################################################
#
# Tests related to performance and memory allocations.
#
############################################################################################

@testset "Aqua.jl" begin
    Aqua.test_all(SatelliteToolboxGeomagneticField; ambiguities=(recursive = false), deps_compat=(check_extras = false))
end

if VERSION >= v"1.12"
    @warn "JET tests skipped on Julia 1.12+ due to JET compatibility limitations"
else
    @testset "JET Testing" begin
        rep = JET.test_package(SatelliteToolboxGeomagneticField; toplevel_logger=nothing, target_modules=(SatelliteToolboxGeomagneticField,))
    end
end

# Skip allocation tests on macOS with Julia 1.12+ due to AllocCheck detecting
# platform-specific runtime calls (jl_get_pgcstack_static) as allocations
if Sys.isapple() && VERSION >= v"1.12"
    @warn "Allocation tests skipped on macOS with Julia 1.12+ due to AllocCheck platform limitations"
else
    @testset "Allocation Check" begin
        @test length(
            check_allocs(
                (date, r, λ, Ω, P, dP) -> begin
                    igrf(date, r, λ, Ω; P = P, dP = dP, verbose = Val(false))
                end,
                (Float64, Float64, Float64, Float64, Matrix{Float64}, Matrix{Float64})
            )
        ) == 0

        @test length(
            check_allocs(
                (date, h, λ, Ω, P, dP) -> begin
                    igrf(date, h, λ, Ω, Val(:geodetic); P = P, dP = dP, verbose = Val(false))
                end,
                (Float64, Float64, Float64, Float64, Matrix{Float64}, Matrix{Float64})
            )
        ) == 0

        @test length(
            check_allocs(
                (date, r, λ, Ω, P, dP) -> begin
                    igrfd(date, r, λ, Ω; P = P, dP = dP, verbose = Val(false))
                end,
                (Float64, Float64, Float64, Float64, Matrix{Float64}, Matrix{Float64})
            )
        ) == 0

        @test length(
            check_allocs(
                (date, h, λ, Ω, P, dP) -> begin
                    igrfd(date, h, λ, Ω, Val(:geodetic); P = P, dP = dP, verbose = Val(false))
                end,
                (Float64, Float64, Float64, Float64, Matrix{Float64}, Matrix{Float64})
            )
        ) == 0

        @test length(
            check_allocs(
                (r_e, year) -> begin
                    geomagnetic_dipole_field(r_e, year)
                end,
                (SVector{3, Float64}, Float64)
            )
        ) == 0

        @test length(
            check_allocs(
                (date, r, λ, Ω, P, dP) -> begin
                    igrf(date, r, λ, Ω; P = P, dP = dP, verbose = Val(false))
                end,
                (Float64, Float64, Float64, Float64, LowerTriangularStorage{RowMajor, Float64}, LowerTriangularStorage{RowMajor, Float64})
            )
        ) == 0

        @test length(
            check_allocs(
                (date, h, λ, Ω, P, dP) -> begin
                    igrf(date, h, λ, Ω, Val(:geodetic); P = P, dP = dP, verbose = Val(false))
                end,
                (Float64, Float64, Float64, Float64, LowerTriangularStorage{RowMajor, Float64}, LowerTriangularStorage{RowMajor, Float64})
            )
        ) == 0

        @test length(
            check_allocs(
                (date, r, λ, Ω, P, dP) -> begin
                    igrfd(date, r, λ, Ω; P = P, dP = dP, verbose = Val(false))
                end,
                (Float64, Float64, Float64, Float64, LowerTriangularStorage{RowMajor, Float64}, LowerTriangularStorage{RowMajor, Float64})
            )
        ) == 0

        @test length(
            check_allocs(
                (date, h, λ, Ω, P, dP) -> begin
                    igrfd(date, h, λ, Ω, Val(:geodetic); P = P, dP = dP, verbose = Val(false))
                end,
                (Float64, Float64, Float64, Float64, LowerTriangularStorage{RowMajor, Float64}, LowerTriangularStorage{RowMajor, Float64})
            )
        ) == 0
    end
end
