## Description #############################################################################
#
# Zygote Extension for the SatelliteToolboxGeomagneticModels.jl package. Needed since igrf
# is mutating and changing that would require a large rework. Instead use ForwardDiff for this function instead.
#
############################################################################################

module SatelliteToolboxGeomagneticModelsZygoteExt

using SatelliteToolboxGeomagneticModels

using Zygote: NoTangent
using Zygote.ChainRulesCore: ChainRulesCore

using ForwardDiff

function ChainRulesCore.rrule(
    ::typeof(GeomagneticModels.igrf),
    date::Number,
    r::AbstractVector{V},
    time::Number;
    max_degree::Int = _IGRF_MAX_DEGREE,
    show_warnings::Bool = true,
    verbosity::Val{verbose} = Val(false),
    P::Union{Nothing, AbstractMatrix} = nothing,
    dP::Union{Nothing, AbstractMatrix} = nothing
) where {V<:Number, verbose}

    y = GeomagneticModels.igrf(
        model,
        r,
        time;
        max_degree = max_degree,
        show_warnings = show_warnings,
        verbosity = verbosity,
        P = P,
        dP = dP,
    )

    function _igrf_pullback(Δ)
        jac = ForwardDiff.jacobian(
            (x) -> GeomagneticModels.igrf(
                date,
                x[1:3],
                x[4];
                max_degree = max_degree,
                show_warnings = show_warnings,
                verbosity = verbosity,
                P = P,
                dP = dP,
            ),
            [r; time]
        )

        vjp = Δ' * jac

        return (NoTangent(), NoTangent(), vjp[1:3], vjp[4], (NoTangent(), NoTangent(), NoTangent()))

    end

    return y, _igrf_pullback

end 

end