## Description #############################################################################
#
# Tests related to performance and memory allocations.
#
############################################################################################

@testset "Aqua.jl" begin
    Aqua.test_all(SatelliteToolboxGeomagneticField; ambiguities = (recursive = false))
end

@testset "JET Testing" begin
    JET.test_package(
        SatelliteToolboxGeomagneticField;
        toplevel_logger = nothing,
        target_modules = (SatelliteToolboxGeomagneticField,),
    )
end

# The allocation-free guarantee of `igrf` and `igrfd` relies on the `@inbounds` annotations
# in the IGRF kernel: the `MVector`s that cache the sines and cosines of the multiples of
# the longitude are stack allocated only if the compiler can remove the bounds-check
# branches, in which the vectors escape to `throw_boundserror`. Hence, the checks are
# meaningless when bounds checking is forced with `--check-bounds=yes`, as `Pkg.test` does
# on Julia < 1.13. CI runs the test suite with `--check-bounds=auto` (see the workflows).
if Base.JLOptions().check_bounds == 1
    @warn "Allocation checks skipped because bounds checking is forced (--check-bounds=yes)"
else
    @testset "Allocation Check" begin
        @test length(
            check_allocs(
                (date, r, λ, Ω, P, dP) -> begin
                    igrf(date, r, λ, Ω; P = P, dP = dP, show_warnings = Val(false))
                end,
                (Float64, Float64, Float64, Float64, Matrix{Float64}, Matrix{Float64}),
            ),
        ) == 0

        @test length(
            check_allocs(
                (date, h, λ, Ω, P, dP) -> begin
                    igrf(
                        date,
                        h,
                        λ,
                        Ω,
                        Val(:geodetic);
                        P = P,
                        dP = dP,
                        show_warnings = Val(false),
                    )
                end,
                (Float64, Float64, Float64, Float64, Matrix{Float64}, Matrix{Float64}),
            ),
        ) == 0

        @test length(
            check_allocs(
                (date, r, λ, Ω, P, dP) -> begin
                    igrfd(date, r, λ, Ω; P = P, dP = dP, show_warnings = Val(false))
                end,
                (Float64, Float64, Float64, Float64, Matrix{Float64}, Matrix{Float64}),
            ),
        ) == 0

        @test length(
            check_allocs(
                (date, h, λ, Ω, P, dP) -> begin
                    igrfd(
                        date,
                        h,
                        λ,
                        Ω,
                        Val(:geodetic);
                        P = P,
                        dP = dP,
                        show_warnings = Val(false),
                    )
                end,
                (Float64, Float64, Float64, Float64, Matrix{Float64}, Matrix{Float64}),
            ),
        ) == 0

        @test length(
            check_allocs(
                (r_e, year) -> begin
                    geomagnetic_dipole_field(r_e, year)
                end,
                (SVector{3, Float64}, Float64),
            ),
        ) == 0

        @test length(
            check_allocs(
                (date, r, λ, Ω, P, dP) -> begin
                    igrf(date, r, λ, Ω; P = P, dP = dP, show_warnings = Val(false))
                end,
                (
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    LowerTriangularStorage{RowMajor, Float64},
                    LowerTriangularStorage{RowMajor, Float64},
                ),
            ),
        ) == 0

        @test length(
            check_allocs(
                (date, h, λ, Ω, P, dP) -> begin
                    igrf(
                        date,
                        h,
                        λ,
                        Ω,
                        Val(:geodetic);
                        P = P,
                        dP = dP,
                        show_warnings = Val(false),
                    )
                end,
                (
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    LowerTriangularStorage{RowMajor, Float64},
                    LowerTriangularStorage{RowMajor, Float64},
                ),
            ),
        ) == 0

        @test length(
            check_allocs(
                (date, r, λ, Ω, P, dP) -> begin
                    igrfd(date, r, λ, Ω; P = P, dP = dP, show_warnings = Val(false))
                end,
                (
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    LowerTriangularStorage{RowMajor, Float64},
                    LowerTriangularStorage{RowMajor, Float64},
                ),
            ),
        ) == 0

        @test length(
            check_allocs(
                (date, h, λ, Ω, P, dP) -> begin
                    igrfd(
                        date,
                        h,
                        λ,
                        Ω,
                        Val(:geodetic);
                        P = P,
                        dP = dP,
                        show_warnings = Val(false),
                    )
                end,
                (
                    Float64,
                    Float64,
                    Float64,
                    Float64,
                    LowerTriangularStorage{RowMajor, Float64},
                    LowerTriangularStorage{RowMajor, Float64},
                ),
            ),
        ) == 0
    end
end
