## Description #############################################################################
#
# Tests for the Zygote/ForwardDiff extension (rrule for igrf geocentric).
#
############################################################################################

@testset "Zygote IGRF Geocentric Gradient" begin
    date = 2020.0
    r    = 6800e3
    λ    = 0.45
    Ω    = -1.34

    function loss_sum(x)
        B = igrf(x[1], x[2], x[3], x[4], Val(:geocentric); show_warnings = false)
        return sum(B)
    end

    g_zyg = Zygote.gradient(loss_sum, [date, r, λ, Ω])[1]
    g_fwd = ForwardDiff.gradient(loss_sum, [date, r, λ, Ω])

    @test g_zyg ≈ g_fwd

    function loss_norm2(x)
        B = igrf(x[1], x[2], x[3], x[4], Val(:geocentric); show_warnings = false)
        return B[1]^2 + B[2]^2 + B[3]^2
    end

    g_zyg = Zygote.gradient(loss_norm2, [date, r, λ, Ω])[1]
    g_fwd = ForwardDiff.gradient(loss_norm2, [date, r, λ, Ω])

    @test g_zyg ≈ g_fwd
end

@testset "Zygote IGRF Geocentric Reduced Degree" begin
    date = 2020.0
    r    = 6800e3
    λ    = 0.30
    Ω    = 0.80

    function loss_reduced(x)
        B = igrf(
            x[1], x[2], x[3], x[4], Val(:geocentric);
            max_degree = 4,
            show_warnings = false,
        )
        return sum(B)
    end

    g_zyg = Zygote.gradient(loss_reduced, [date, r, λ, Ω])[1]
    g_fwd = ForwardDiff.gradient(loss_reduced, [date, r, λ, Ω])

    @test g_zyg ≈ g_fwd
end

@testset "Zygote IGRF Geocentric Different Locations" begin
    test_cases = [
        (2015.0, 7000e3,  0.78, -0.50),
        (2005.5, 6500e3, -0.30,  2.10),
        (2025.0, 6900e3,  1.20, -2.80),
    ]

    for (date, r, λ, Ω) in test_cases
        function loss(x)
            B = igrf(x[1], x[2], x[3], x[4], Val(:geocentric); show_warnings = false)
            return sum(B)
        end

        x = [date, r, λ, Ω]
        g_zyg = Zygote.gradient(loss, x)[1]
        g_fwd = ForwardDiff.gradient(loss, x)

        @test g_zyg ≈ g_fwd
    end
end
