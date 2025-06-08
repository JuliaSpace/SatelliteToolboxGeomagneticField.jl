## Description #############################################################################
#
# Tests related to the IGRF model.
#
############################################################################################

# Load test files.
#
# The test values were created using the original `igrf14syn` function. The source code of
# the program that generated those files can be seen at:
#
#   https://github.com/JuliaSpace/IGRF_Test

const igrf14_geocentric_test = readdlm("./IGRF14_test_geocentric.txt")
const igrf14_geodetic_test   = readdlm("./IGRF14_test_geodetic.txt")

# == File: ./src/igrf/igrf.jl ==============================================================

# -- Function: igrf ------------------------------------------------------------------------

@testset "Function: igrf" begin
    # Auxiliary variables to use the version without allocations.
    P  = Matrix{Float64}(undef, 14, 14)
    dP = similar(P)

    # Testing the geocentric part of the algorithm.
    for i in 1:size(igrf14_geocentric_test, 1)
        date  = igrf14_geocentric_test[i, 1]
        r     = igrf14_geocentric_test[i, 2]
        colat = igrf14_geocentric_test[i, 3]
        elong = igrf14_geocentric_test[i, 4]
        xt    = igrf14_geocentric_test[i, 5]
        yt    = igrf14_geocentric_test[i, 6]
        zt    = igrf14_geocentric_test[i, 7]
        ft    = igrf14_geocentric_test[i, 8]

        # Call IGRF with the same inputs as those in the test.
        if elong > 180
            elong = elong - 360
        end

        if date <= 2030
            Ba = igrf(date, 1000r, deg2rad(90 - colat), deg2rad(elong))
            Bn = igrf(date, 1000r, deg2rad(90 - colat), deg2rad(elong); P = P, dP = dP)
        else
            Ba = @test_logs (
                :warn,
                "The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than 2030."
            ) igrf(date, 1000r, deg2rad(90 - colat), deg2rad(elong))

            Bn = @test_logs (
                :warn,
                "The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than 2030."
            ) igrf(date, 1000r, deg2rad(90 - colat), deg2rad(elong); P = P, dP = dP)
        end

        x, y, z = Ba[:]
        f       = norm(Ba)

        # Test the values.
        @test x   ≈ xt atol = 3e-1
        @test y   ≈ yt atol = 3e-1
        @test z   ≈ zt atol = 3e-1
        @test f   ≈ ft atol = 3e-1
        @test Ba == Bn
    end

    # Testing the geodetic part of the algorithm.
    for i in 1:size(igrf14_geocentric_test, 1)
        date  = igrf14_geodetic_test[i, 1]
        h     = igrf14_geodetic_test[i, 2]
        colat = igrf14_geodetic_test[i, 3]
        elong = igrf14_geodetic_test[i, 4]
        xt    = igrf14_geodetic_test[i, 5]
        yt    = igrf14_geodetic_test[i, 6]
        zt    = igrf14_geodetic_test[i, 7]
        ft    = igrf14_geodetic_test[i, 8]

        # Call IGRF with the same inputs as those in the test.
        if elong > 180
            elong = elong - 360
        end

        if date <= 2030
            Ba = igrf(date, 1000h, deg2rad(90 - colat), deg2rad(elong), Val(:geodetic))
            Bn = igrf(date, 1000h, deg2rad(90 - colat), deg2rad(elong), Val(:geodetic); P = P, dP = dP)
        else
            Ba = @test_logs (
                :warn,
                "The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than 2030."
            ) igrf(date, 1000h, deg2rad(90 - colat), deg2rad(elong), Val(:geodetic))

            Bn = @test_logs (
                :warn,
                "The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than 2030."
            ) igrf(date, 1000h, deg2rad(90 - colat), deg2rad(elong), Val(:geodetic); P = P, dP = dP)
        end

        x, y, z = Ba[:]
        f       = norm(Ba)

        # Test the values.
        @test x   ≈ xt atol = 3e-1
        @test y   ≈ yt atol = 3e-1
        @test z   ≈ zt atol = 3e-1
        @test f   ≈ ft atol = 3e-1
        @test Ba == Bn
    end

    # == Reduced Degree ====================================================================

    # -- Geocentric ------------------------------------------------------------------------

    # If `max_degree` is equal to or lower than 0, the results must be the same as if
    # `max_degree` is 1.
    Bref = igrf(2020.4452, 6515e3, 0.45, -1.34; max_degree =  1)
    B1   = igrf(2020.4452, 6515e3, 0.45, -1.34; max_degree =  0)
    B2   = igrf(2020.4452, 6515e3, 0.45, -1.34; max_degree = -1)

    @test Bref == B1
    @test Bref == B2

    # If `max_degree` is higher than the number of available coefficients, it must be
    # clamped.
    Bref = igrf(2020.4452, 6515e3, 0.45, -1.34)
    B1   = igrf(2020.4452, 6515e3, 0.45, -1.34; max_degree = 13)
    B2   = igrf(2020.4452, 6515e3, 0.45, -1.34; max_degree = 22)

    @test Bref == B1
    @test Bref == B2

    # The difference between the field computed using 13 or 12 coefficients must be small.
    B1 = igrf(2020.4452, 6515e3, 0.45, -1.34)
    B2 = igrf(2020.4452, 6515e3, 0.45, -1.34; max_degree = 12)

    @test acosd(dot(B1 / norm(B1), B2 / norm(B2))) < 0.01

    # Testing the version in which we provide the matrices to the associated Legendre
    # functions.
    Pred  = zeros(5, 5)
    dPred = zeros(5, 5)

    Bref = igrf(2020.4452, 6515e3, 0.45, -1.34; max_degree = 4)
    B1   = igrf(2020.4452, 6515e3, 0.45, -1.34; P = Pred, dP = dPred, max_degree = 4)

    @test B1 == Bref

    # -- Geodetic --------------------------------------------------------------------------

    # If `max_degree` is equal or lower than 0, the results must be the same as if
    # `max_degree` is 1.
    Bref = igrf(2020.4452, 752e3, 0.45, -1.34, Val(:geodetic); max_degree =  1)
    B1   = igrf(2020.4452, 752e3, 0.45, -1.34, Val(:geodetic); max_degree =  0)
    B2   = igrf(2020.4452, 752e3, 0.45, -1.34, Val(:geodetic); max_degree = -1)

    @test Bref == B1
    @test Bref == B2

    # If `max_degree` is higher than the number of available coefficients, it must be
    # clamped.
    Bref = igrf(2020.4452, 752e3, 0.45, -1.34, Val(:geodetic))
    B1   = igrf(2020.4452, 752e3, 0.45, -1.34, Val(:geodetic); max_degree = 13)
    B2   = igrf(2020.4452, 752e3, 0.45, -1.34, Val(:geodetic); max_degree = 22)

    @test Bref == B1
    @test Bref == B2

    # The difference between the field computed using 13 or 12 coefficients must be small.
    B1 = igrf(2020.4452, 752e3, 0.45, -1.34, Val(:geodetic))
    B2 = igrf(2020.4452, 752e3, 0.45, -1.34, Val(:geodetic); max_degree = 12)

    @test acosd(dot(B1 / norm(B1), B2 / norm(B2))) < 0.01

    # Testing the version in which we provide the matrices to the associated Legendre
    # functions.
    Pred  = zeros(5, 5)
    dPred = zeros(5, 5)

    Bref = igrf(2020.4452, 752e3, 0.45, -1.34, Val(:geodetic); max_degree = 4)
    B1   = igrf(2020.4452, 752e3, 0.45, -1.34, Val(:geodetic); P = Pred, dP = dPred, max_degree = 4)

    @test B1 == Bref

    # == Float32 ===========================================================================

    Bf32 = igrf(2020.4452, 6515f3, 0.45f0, -1.34f0)
    @test eltype(Bf32) === Float32

    Bf32 = igrf(2020.4452, 6515f3, 0, -1)
    @test eltype(Bf32) === Float32
end

@testset "Function: igrf [Issues]" begin
    # Calculation close to the geographic pole.
    B = igrf(2019, 7150e3, π / 2 - 1e-15, 0.55)
    @test B[1] ≈   908.1899663358307
    @test B[2] ≈   173.01468080386584
    @test B[3] ≈ 41139.84620809845

    B = igrf(2019, 7150e3, π / 2, 0.55)
    @test B[1] ≈   908.1899663358307
    @test B[2] ≈   173.01468080386584
    @test B[3] ≈ 41139.84620809845
end

@testset "Function: igrf [ERRORS]" begin
    P₀  = zeros(10, 10)
    dP₀ = zeros(10, 10)
    P₁  = zeros(15, 15)
    dP₁ = zeros(15, 15)
    R0  = 6378.137e3

    @test_throws ArgumentError igrf(2020, 140e3, π / 4, π / 2, Val(:geodetic); P = P₀, dP = dP₀)
    @test_throws ArgumentError igrf(2020, 140e3, π / 4, π / 2, Val(:geodetic); P = P₀, dP = dP₁)
    @test_throws ArgumentError igrf(2020, 140e3, π / 4, π / 2, Val(:geodetic); P = P₁, dP = dP₀)
    @test_nowarn igrf(2020, 140e3, π / 4, π / 2, Val(:geodetic); P = P₁, dP = dP₁)

    @test_throws ArgumentError igrf(2020, R0 + 140e3, π / 4, π / 2, Val(:geocentric); P = P₀, dP = dP₀)
    @test_throws ArgumentError igrf(2020, R0 + 140e3, π / 4, π / 2, Val(:geocentric); P = P₀, dP = dP₁)
    @test_throws ArgumentError igrf(2020, R0 + 140e3, π / 4, π / 2, Val(:geocentric); P = P₁, dP = dP₀)
    @test_nowarn igrf(2020, R0 + 140e3, π / 4, π / 2, Val(:geocentric); P = P₁, dP = dP₁)

    @test_throws ArgumentError igrf(1899.9, R0 + 140e3, π / 4, π / 2, Val(:geocentric))
    @test_throws ArgumentError igrf(2035.1, R0 + 140e3, π / 4, π / 2, Val(:geocentric))

    @test_throws ArgumentError igrf(2020, R0 + 140e3, -π / 2 - 0.01, π / 2, Val(:geocentric))
    @test_throws ArgumentError igrf(2020, R0 + 140e3, +π / 2 + 0.01, π / 2, Val(:geocentric))
    @test_throws ArgumentError igrf(2020, R0 + 140e3, π / 4, -π - 0.01, Val(:geocentric))
    @test_throws ArgumentError igrf(2020, R0 + 140e3, π / 4, +π + 0.01, Val(:geocentric))
end

# -- Function: igrfd -----------------------------------------------------------------------

@testset "Function: igrfd" begin
    # Auxiliary variables to use the version without allocations.
    P  = Matrix{Float64}(undef, 14, 14)
    dP = similar(P)

    # Testing the geocentric part of the algorithm.
    for i in 1:size(igrf14_geocentric_test, 1)
        date  = igrf14_geocentric_test[i, 1]
        r     = igrf14_geocentric_test[i, 2]
        colat = igrf14_geocentric_test[i, 3]
        elong = igrf14_geocentric_test[i, 4]
        xt    = igrf14_geocentric_test[i, 5]
        yt    = igrf14_geocentric_test[i, 6]
        zt    = igrf14_geocentric_test[i, 7]
        ft    = igrf14_geocentric_test[i, 8]

        # Call IGRF with the same inputs as those in the test.
        if elong > 180
            elong = elong - 360
        end

        if date <= 2030
            Ba = igrfd(date, 1000r, 90 - colat, elong)
            Bn = igrfd(date, 1000r, 90 - colat, elong; P = P, dP = dP)
        else
            Ba = @test_logs (
                :warn,
                "The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than 2030."
            ) igrfd(date, 1000r, 90 - colat, elong)

            Bn = @test_logs (
                :warn,
                "The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than 2030."
            ) igrfd(date, 1000r, 90 - colat, elong; P = P, dP = dP)
        end

        x, y, z = Ba[:]
        f       = norm(Ba)

        # Test the values.
        @test x   ≈ xt atol = 3e-1
        @test y   ≈ yt atol = 3e-1
        @test z   ≈ zt atol = 3e-1
        @test f   ≈ ft atol = 3e-1
        @test Ba == Bn
    end

    # Testing the geodetic part of the algorithm.
    for i in 1:size(igrf14_geocentric_test, 1)
        date  = igrf14_geodetic_test[i, 1]
        h     = igrf14_geodetic_test[i, 2]
        colat = igrf14_geodetic_test[i, 3]
        elong = igrf14_geodetic_test[i, 4]
        xt    = igrf14_geodetic_test[i, 5]
        yt    = igrf14_geodetic_test[i, 6]
        zt    = igrf14_geodetic_test[i, 7]
        ft    = igrf14_geodetic_test[i, 8]

        # Call IGRF with the same inputs as those in the test.
        if elong > 180
            elong = elong - 360
        end

        if date <= 2030
            Ba = igrfd(date, 1000h, 90 - colat, elong, Val(:geodetic))
            Bn = igrfd(date, 1000h, 90 - colat, elong, Val(:geodetic); P = P, dP = dP)
        else
            Ba = @test_logs (
                :warn,
                "The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than 2030."
            ) igrfd(date, 1000h, 90 - colat, elong, Val(:geodetic))

            Bn = @test_logs (
                :warn,
                "The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than 2030."
            ) igrfd(date, 1000h, 90 - colat, elong, Val(:geodetic); P = P, dP = dP)
        end

        x, y, z = Ba[:]
        f       = norm(Ba)

        # Test the values.
        @test x   ≈ xt atol = 3e-1
        @test y   ≈ yt atol = 3e-1
        @test z   ≈ zt atol = 3e-1
        @test f   ≈ ft atol = 3e-1
        @test Ba == Bn
    end

    # == Reduced Degree ====================================================================

    # -- Geocentric ------------------------------------------------------------------------

    # If `max_degree` is equal to or lower than 0, the results must be the same as if
    # `max_degree` is 1.
    Bref = igrfd(2020.4452, 6515e3, -19, 106; max_degree =  1)
    B1   = igrfd(2020.4452, 6515e3, -19, 106; max_degree =  0)
    B2   = igrfd(2020.4452, 6515e3, -19, 106; max_degree = -1)

    @test Bref == B1
    @test Bref == B2

    # If `max_degree` is higher than the number of available coefficients, it must be
    # clamped.
    Bref = igrfd(2020.4452, 6515e3, -19, 106)
    B1   = igrfd(2020.4452, 6515e3, -19, 106; max_degree = 13)
    B2   = igrfd(2020.4452, 6515e3, -19, 106; max_degree = 22)

    @test Bref == B1
    @test Bref == B2

    # The difference between the field computed using 13 or 12 coefficients must be small.
    B1 = igrfd(2020.4452, 6515e3, -19, 106)
    B2 = igrfd(2020.4452, 6515e3, -19, 106; max_degree = 12)

    @test acosd(dot(B1 / norm(B1), B2 / norm(B2))) < 0.01

    # Testing the version in which we provide the matrices to the associated Legendre
    # functions.
    Pred  = zeros(5, 5)
    dPred = zeros(5, 5)

    Bref = igrfd(2020.4452, 6515e3, -19, 106; max_degree = 4)
    B1   = igrfd(2020.4452, 6515e3, -19, 106; P = Pred, dP = dPred, max_degree = 4)

    @test B1 == Bref

    # -- Geodetic --------------------------------------------------------------------------

    # If `max_degree` is equal or lower than 0, the results must be the same as if
    # `max_degree` is 1.
    Bref = igrfd(2020.4452, 752e3, -19, 106, Val(:geodetic); max_degree =  1)
    B1   = igrfd(2020.4452, 752e3, -19, 106, Val(:geodetic); max_degree =  0)
    B2   = igrfd(2020.4452, 752e3, -19, 106, Val(:geodetic); max_degree = -1)

    @test Bref == B1
    @test Bref == B2

    # If `max_degree` is higher than the number of available coefficients, it must be
    # clamped.
    Bref = igrfd(2020.4452, 752e3, -19, 106, Val(:geodetic))
    B1   = igrfd(2020.4452, 752e3, -19, 106, Val(:geodetic); max_degree = 13)
    B2   = igrfd(2020.4452, 752e3, -19, 106, Val(:geodetic); max_degree = 22)

    @test Bref == B1
    @test Bref == B2

    # The difference between the field computed using 13 or 12 coefficients must be small.
    B1 = igrfd(2020.4452, 752e3, -19, 106, Val(:geodetic))
    B2 = igrfd(2020.4452, 752e3, -19, 106, Val(:geodetic); max_degree = 12)

    @test acosd(dot(B1 / norm(B1), B2 / norm(B2))) < 0.01

    # Testing the version in which we provide the matrices to the associated Legendre
    # functions.
    Pred  = zeros(5, 5)
    dPred = zeros(5, 5)

    Bref = igrfd(2020.4452, 752e3, -19, 106, Val(:geodetic); max_degree = 4)
    B1   = igrfd(2020.4452, 752e3, -19, 106, Val(:geodetic); P = Pred, dP = dPred, max_degree = 4)

    @test B1 == Bref

    # == Float32 ===========================================================================

    Bf32 = igrfd(2020.4452, 6515f3, -25.0f0, -45.0f0)
    @test eltype(Bf32) === Float32

    Bf32 = igrfd(2020.4452, 6515f3, -25, -45)
    @test eltype(Bf32) === Float32
end

@testset "Function: igrfd [ERRORS]" begin
    P₀  = zeros(10, 10)
    dP₀ = zeros(10, 10)
    P₁  = zeros(15, 15)
    dP₁ = zeros(15, 15)
    R0  = 6378.137e3

    @test_throws ArgumentError igrfd(2020, 140e3, 45, 90, Val(:geodetic); P = P₀, dP = dP₀)
    @test_throws ArgumentError igrfd(2020, 140e3, 45, 90, Val(:geodetic); P = P₀, dP = dP₁)
    @test_throws ArgumentError igrfd(2020, 140e3, 45, 90, Val(:geodetic); P = P₁, dP = dP₀)
    @test_nowarn igrfd(2020, 140e3, 45, 90, Val(:geodetic); P = P₁, dP = dP₁)

    @test_throws ArgumentError igrfd(2020, R0 + 140e3, 45, 90, Val(:geocentric); P = P₀, dP = dP₀)
    @test_throws ArgumentError igrfd(2020, R0 + 140e3, 45, 90, Val(:geocentric); P = P₀, dP = dP₁)
    @test_throws ArgumentError igrfd(2020, R0 + 140e3, 45, 90, Val(:geocentric); P = P₁, dP = dP₀)
    @test_nowarn igrfd(2020, R0 + 140e3, 45, 90, Val(:geocentric); P = P₁, dP = dP₁)

    @test_throws ArgumentError igrfd(2020, R0 + 140e3, -91,   90, Val(:geocentric))
    @test_throws ArgumentError igrfd(2020, R0 + 140e3, +91,   90, Val(:geocentric))
    @test_throws ArgumentError igrfd(2020, R0 + 140e3, +45,  181, Val(:geocentric))
    @test_throws ArgumentError igrfd(2020, R0 + 140e3, +45, -181, Val(:geocentric))
end
