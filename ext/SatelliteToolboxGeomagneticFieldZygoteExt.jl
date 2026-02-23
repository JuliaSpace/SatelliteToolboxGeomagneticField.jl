## Description #############################################################################
#
# Zygote Extension for the SatelliteToolboxGeomagneticField.jl package. Needed since igrf
# is mutating and changing that would require a large rework. Instead use ForwardDiff for this function instead.
#
############################################################################################

module SatelliteToolboxGeomagneticFieldZygoteExt

using SatelliteToolboxGeomagneticField
using SatelliteToolboxGeomagneticField: _IGRF_MAX_DEGREE

using Zygote: NoTangent
using Zygote.ChainRulesCore: ChainRulesCore

using ForwardDiff

function ChainRulesCore.rrule(
    ::typeof(igrf),
    date::Number,
    r::Number,
    λ_gc::Number,
    Ω::Number,
    ::Val{:geocentric};
    max_degree::Int = _IGRF_MAX_DEGREE,
    show_warnings::Bool = true,
    verbosity::Val{verbose} = Val(false),
    P::Union{Nothing, AbstractMatrix} = nothing,
    dP::Union{Nothing, AbstractMatrix} = nothing
) where {verbose}

    y = igrf(
        date,
        r,
        λ_gc,
        Ω,
        Val(:geocentric);
        max_degree = max_degree,
        show_warnings = show_warnings,
        verbosity = verbosity,
        P = P,
        dP = dP,
    )

    function _igrf_pullback(Δ)
        jac = ForwardDiff.jacobian(
            (x) -> igrf(
                x[1],
                x[2],
                x[3],
                x[4],
                Val(:geocentric);
                max_degree = max_degree,
                show_warnings = show_warnings,
                verbosity = verbosity,
            ),
            [date; r; λ_gc; Ω]
        )

        vjp = Δ' * jac

        return (NoTangent(), vjp[1], vjp[2], vjp[3], vjp[4], NoTangent())
    end

    return y, _igrf_pullback

end

end