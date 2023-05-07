# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   International Geomagnetic Field Model.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
#   [2] https://www.ngdc.noaa.gov/IAGA/vmod/igrf12.f
#   [3] https://www.mathworks.com/matlabcentral/fileexchange/34388-international-geomagnetic-reference-field--igrf--model
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export igrf, igrfd

############################################################################################
#                                        Functions
############################################################################################

"""
    igrfd(date::Number, <r, h>::T1, λ::T2, Ω::T3[, R]; kwargs...) where {T1<:Number, T2<:Number, T3<:Number} -> SVector{3, T}

**IGRF Model**

*Current version: v13*

Compute the geomagnetic field vector [nT] at the date `date` [Year A.D.] and position (`r`
or `h`, `λ`, `Ω`).

The position representation is defined by `R`. If `R` is `Val(:geocentric)`, the input must
be **geocentric** coordinates:

1. Distance from the Earth center `r` [m];
2. Geocentric latitude `λ` ∈ (-90°, +90°); and
3. Geocentric longitude `Ω` ∈ (-180°, +180°).

If `R` is `Val(:geodetic)`, the input must be **geodetic** coordinates:

1 Altitude above the reference ellipsoid `h` (WGS-84) [m];
2. Geodetic latitude `λ` ∈ (-90°, +90°); and
3. Geodetic longitude `Ω` ∈ (-180°, +180°).

If `R` is omitted, it defaults to `Val(:geocentric)`.

!!! warning
    We must have `1900 <= date <= 2030`. A warning message is printed for dates greater than
    2025 since the output is not reliable anymore. This message can be suppressed by setting
    the keyword `show_warnings` to `false`.

!!! info
    The output vector will be represented in the same reference system selected by the
    parameter `R` (geocentric or geodetic). The Y-axis of the output reference system always
    points East. In case of **geocentric coordinates**, the Z-axis points toward the center
    of Earth and the X-axis completes a right-handed coordinate system. In case of
    **geodetic coordinates**, the X-axis is tangent to the ellipsoid at the selected
    location and points toward North, whereas the Z-axis completes a right-hand coordinate
    system.

# Keywords

- `max_degree::Int`: Maximum degree used in the spherical harmonics when computing the
    geomagnetic field. If it is higher than the available number of coefficients in the IGRF
    matrices, it will be clamped. If it is equal of lower than 0, it will be set to 1.
    (**Default** = 13)
- `show_warnings::Bool`: Show warnings about the data (**Default** = `true`).
- `P::Union{Nothing, AbstractMatrix}`: An optional matrix that must contain at least
    `max_degree + 1 × max_degree + 1` real numbers that will be used to store the Legendre
    coefficients, reducing the allocations. If it is `nothing`, the matrix will be created
    when calling the function.
- `dP::Union{Nothing, AbstractMatrix}`: An optional matrix that must contain at least
    `max_degree + 1 × max_degree + 1` real numbers that will be used to store the Legendre
    derivative coefficients, reducing the allocations. If it is `nothing`, the matrix will
    be created when calling the function.

# Returns

- `SVector{3, T}`: Geomagnetic field vector [nT] at the desired location represented in the
    same input reference (geocentric or geodetic).

!!! info
    The output type `T` is obtained by promoting `T1`, `T2`, and `T3` to a float.
"""
function igrfd(
    date::Number,
    r::Number,
    λ_gc::Number,
    Ω::Number;
    max_degree::Int = _IGRF_MAX_DEGREE,
    show_warnings::Bool = true,
    P::Union{Nothing, AbstractMatrix} = nothing,
    dP::Union{Nothing, AbstractMatrix} = nothing
)
    return igrfd(
        date,
        r,
        λ_gc,
        Ω,
        Val(:geocentric);
        max_degree = max_degree,
        show_warnings = show_warnings,
        P = P,
        dP = dP
    )
end

function igrfd(
    date::Number,
    r::T1,
    λ_gc::T2,
    Ω::T3,
    ::Val{:geocentric};
    max_degree::Int = _IGRF_MAX_DEGREE,
    show_warnings::Bool = true,
    P::Union{Nothing, AbstractMatrix} = nothing,
    dP::Union{Nothing, AbstractMatrix} = nothing
) where {T1<:Number, T2<:Number, T3<:Number}

    T = promote_type(T1, T2, T3) |> float

    # Check if the latitude and longitude are valid.
    if (λ_gc < -90) || (λ_gc > 90)
        throw(ArgumentError("The latitude must be between -90° and +90° rad."))
    end

    if (Ω < -180) || (Ω > 180)
        throw(ArgumentError("The longitude must be between -180° and +180° rad."))
    end

    return igrf(
        date,
        T(r),
        deg2rad(T(λ_gc)),
        deg2rad(T(Ω)),
        Val(:geocentric);
        max_degree = max_degree,
        show_warnings = show_warnings,
        P = P,
        dP = dP
    )
end

function igrfd(
    date::Number,
    h::T1,
    λ::T2,
    Ω::T3,
    ::Val{:geodetic};
    max_degree::Int = _IGRF_MAX_DEGREE,
    show_warnings::Bool = true,
    P::Union{Nothing, AbstractMatrix} = nothing,
    dP::Union{Nothing, AbstractMatrix} = nothing
) where {T1<:Number, T2<:Number, T3<:Number}

    T = promote_type(T1, T2, T3) |> float

    # Check if the latitude and longitude are valid.
    if (λ < -90) || (λ > 90)
        throw(ArgumentError("The latitude must be between -90° and +90° rad."))
    end

    if (Ω < -180) || (Ω > 180)
        throw(ArgumentError("The longitude must be between -180° and +180° rad."))
    end

    return igrf(
        date,
        T(h),
        deg2rad(T(λ)),
        deg2rad(T(Ω)),
        Val(:geodetic);
        max_degree = max_degree,
        show_warnings = show_warnings,
        P = P,
        dP = dP
    )
end

"""
    igrf(date::Number, <r, h>::T1, λ::T2, Ω::T3[, R]; kwargs...) where {T1<:Number, T2<:Number, T3<:Number} -> SVector{3, T}

**IGRF Model**

*Current version: v13*

Compute the geomagnetic field vector [nT] at the date `date` [Year A.D.] and position (`r`
or `h`, `λ`, `Ω`).

The position representation is defined by `R`. If `R` is `Val(:geocentric)`, the input must
be **geocentric** coordinates:

1. Distance from the Earth center `r` [m];
2. Geocentric latitude `λ` ∈ (-π / 2, +π / 2) [rad]; and
3. Geocentric longitude `Ω` ∈ (-π, +π) [rad].

If `R` is `Val(:geodetic)`, the input must be **geodetic** coordinates:

1. Altitude above the reference ellipsoid `h` (WGS-84) \\[m];
2. Geodetic latitude `λ` ∈ (-π/2, +π/2) [rad]; and
3. Geodetic longitude `Ω` ∈ (-π, +π) [rad].

If `R` is omitted, it defaults to `Val(:geocentric)`.

!!! warning
    We must have `1900 <= date <= 2030`. A warning message is printed for dates greater than
    2025 since the output is not reliable anymore. This message can be suppressed by setting
    the keyword `show_warnings` to `false`.

!!! info
    The output vector will be represented in the same reference system selected by the
    parameter `R` (geocentric or geodetic). The Y-axis of the output reference system always
    points East. In case of **geocentric coordinates**, the Z-axis points toward the center
    of Earth and the X-axis completes a right-handed coordinate system. In case of
    **geodetic coordinates**, the X-axis is tangent to the ellipsoid at the selected
    location and points toward North, whereas the Z-axis completes a right-hand coordinate
    system.

# Keywords

- `max_degree::Int`: Maximum degree used in the spherical harmonics when computing the
    geomagnetic field. If it is higher than the available number of coefficients in the IGRF
    matrices, it will be clamped. If it is equal of lower than 0, it will be set to 1.
    (**Default** = 13)
- `show_warnings::Bool`: Show warnings about the data (**Default** = `true`).
- `P::Union{Nothing, AbstractMatrix}`: An optional matrix that must contain at least
    `max_degree + 1 × max_degree + 1` real numbers that will be used to store the Legendre
    coefficients, reducing the allocations. If it is `nothing`, the matrix will be created
    when calling the function.
- `dP::Union{Nothing, AbstractMatrix}`: An optional matrix that must contain at least
    `max_degree + 1 × max_degree + 1` real numbers that will be used to store the Legendre
    derivative coefficients, reducing the allocations. If it is `nothing`, the matrix will
    be created when calling the function.

# Returns

- `SVector{3, T}`: Geomagnetic field vector [nT] at the desired location represented in the
    same input reference (geocentric or geodetic).

!!! info
    The output type `T` is obtained by promoting `T1`, `T2`, and `T3` to a float.
"""
function igrf(
    date::Number,
    r::Number,
    λ_gc::Number,
    Ω::Number;
    max_degree::Int = _IGRF_MAX_DEGREE,
    show_warnings::Bool = true,
    P::Union{Nothing, AbstractMatrix} = nothing,
    dP::Union{Nothing, AbstractMatrix} = nothing
)
    return igrf(
        date,
        r,
        λ_gc,
        Ω,
        Val(:geocentric);
        max_degree = max_degree,
        show_warnings = show_warnings,
        P = P,
        dP = dP
    )
end

function igrf(
    date::Number,
    r::T1,
    λ::T2,
    Ω::T3,
    ::Val{:geocentric};
    max_degree::Int = _IGRF_MAX_DEGREE,
    show_warnings::Bool = true,
    P::Union{Nothing, AbstractMatrix} = nothing,
    dP::Union{Nothing, AbstractMatrix} = nothing
) where {T1<:Number, T2<:Number, T3<:Number}

    T = promote_type(T1, T2, T3) |> float

    # Input Verification
    # ======================================================================================

    # Check the data, since this model is valid for years between 1900 and `max_year`.
    if (date < 1900) || (date > _IGRF_LAST_YEAR)
        throw(ArgumentError(
            "This IGRF version will not work for years outside the interval [1900, $_IGRF_LAST_YEAR)."
        ))
    end

    # Check if the latitude and longitude are valid.
    if (λ < -π/2) || (λ > π/2)
        throw(ArgumentError("The latitude must be between -π / 2 and +π / 2 rad."))
    end

    if (Ω < -π) || (Ω > π)
        throw(ArgumentError("The longitude must be between -π and +π rad."))
    end

    # Warn the user that for dates after the year `rel_year` the accuracy maybe reduced.
    if show_warnings && (date > _IGRF_RELIABLE_YEAR)
        @warn("The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than $_IGRF_RELIABLE_YEAR.")
    end

    # If the `max_degree` is equal or lower than 0, we must clamp it to 1.
    max_degree = max(max_degree, 1)

    # Input Variables Conversion
    # ======================================================================================

    # Convert latitude / longitude to co-latitude and east-longitude.
    θ = T(π / 2) - T(λ)
    ϕ = (Ω >= 0) ? T(Ω) : T(2π) + T(Ω)

    # The input variable `r` is in [m], but all the algorithm requires it to be in [km].
    r_km = T(r) / 1000

    # Preliminary Setup
    # ======================================================================================

    # Compute the epoch that will be used to obtain the coefficients. This is necessary
    # because the IGRF provides coefficients every 5 years. Between two epochs, those
    # coefficients must be interpolated.
    idx   = floor(Int, clamp((date - 1900) / 5 + 1, 0, (_IGRF_RELIABLE_YEAR - 1900) / 5))
    epoch = 1900 + (idx - 1) * 5

    # We must jump the first two columns that are reserved for the degree and order.
    idx += 2

    # Compute the fraction of time from the epoch of the coefficient selected by `idx`.
    Δt = T(date - epoch)

    # Compute the maximum spherical harmonic degree for the selected date.
    n_max = (epoch < 1995) ? 10 : 13

    # Check if the user wants a lower degree.
    n_max = min(max_degree, n_max)

    # Check if the matrices related to Legendre must be computed.
    if isnothing(P)
        P = Matrix{T}(undef, n_max + 1, n_max + 1)

    else
        # If the user passed a matrix, we must check if there are enough space to store the
        # coefficients.
        rows, cols = size(P)

        if (rows < n_max + 1) || (cols < n_max + 1)
            throw(ArgumentError("Matrix `P` must have at least $(n_max + 1) rows and columns."))
        end
    end

    if isnothing(dP)
        dP = Matrix{T}(undef, n_max + 1, n_max + 1)

    else
        # If the user passed a matrix, we must check if there are enough space to store the
        # coefficients.
        rows, cols = size(dP)

        if (rows < n_max + 1) || (cols < n_max + 1)
            throw(ArgumentError("Matrix `P` must have at least $(n_max + 1) rows and columns."))
        end
    end

    # Geomagnetic potential gradient
    # ======================================================================================

    dVr, dVϕ, dVθ = _igrf_geomagnetic_potential_gradient(
        n_max,
        idx,
        r_km,
        θ,
        ϕ,
        Δt,
        date >= _IGRF_LAST_YEAR_WITH_MEASUREMENTS,
        P,
        dP,
    )

    # Compute the Geomagnetic field vector in the geocentric reference frame
    # ======================================================================================

    x = +dVθ / r_km
    y = (θ == 0) ? -dVϕ / r_km : -dVϕ / (r_km * sin(θ))
    z = dVr

    B_gc = SVector{3, T}(x, y, z)

    return B_gc
end

function igrf(
    date::Number,
    h::Number,
    λ::Number,
    Ω::Number,
    ::Val{:geodetic};
    max_degree::Int = _IGRF_MAX_DEGREE,
    show_warnings::Bool = true,
    P::Union{Nothing, AbstractMatrix} = nothing,
    dP::Union{Nothing, AbstractMatrix} = nothing
)

    # TODO: This method has a small error (≈ 0.01 nT) compared with the `igrf12syn`.
    # However, the result is exactly the same as the MATLAB function in [3]. Hence, this
    # does not seem to be an error in the conversion from geodetic to geocentric
    # coordinates. This is probably caused by a numerical error. Further verification is
    # necessary.

    # Convert the geodetic coordinates to geocentric coordinates.
    λ_gc, r = geodetic_to_geocentric(λ, h)

    # Compute the geomagnetic field in geocentric coordinates.
    B_gc = igrf(
        date,
        r,
        λ_gc,
        Ω,
        Val(:geocentric);
        max_degree = max_degree,
        show_warnings = show_warnings,
        P = P,
        dP = dP
    )

    # Convert to geodetic coordinates.
    D_gd_gc = angle_to_dcm(λ_gc - λ, :Y)
    B_gd    = D_gd_gc * B_gc

    return B_gd
end

############################################################################################
#                                    Private Functions
############################################################################################

"""
    _igrf_geomagnetic_potential_gradient(n_max::Int, idx::Int, r_km::T, θ::T, ϕ::T, Δt::T, extrapolate::Bool, P::AbstractMatrix, dP::AbstractMatrix) where T<:Number -> NTuple{3, T}

Compute the geomagnetic potential gradient.

# Arguments

- `n_max::Int`: Maximum degree when computing the potential in the spherical harmonics.
- `idx::Int: Index related to the desired epoch in the matrices `_IGRF_G` and `_IGRF_H`.
- `r_km::T`: Position to compute the gradient from the Earth's center [km].
- `θ::T`: Geocentric co-latitude [rad] ∈ [0, π].
- `ϕ::T`: East-longitude [rad] ∈ [0, 2π].
- `Δt::T`: Elapsed time from the epoch related to the index `idx` [year].
- `extrapolate::Bool`: If `true`, the desired epoch is after the last year with
    measurements (`_IGRF_LAST_YEAR_WITH_MEASUREMENTS`). Hence, we must use the coefficient
    time-derivative in the last column of the matrices `_IGRF_G` and `_IGRF_H`.
- `P::AbstractMatrix`: An auxiliary matrix to compute the values of the Legendre associated
    functions. It must have a dimension equal to or greater than `n_max + 1 × n_max + 1`.
- `dP::AbstractMatrix`: An auxiliary matrix to compute the derivatives of the Legendre
    associated functions. It must have a dimension equal to or greater than
    `n_max + 1 × n_max + 1`.

!!! warning
    This is a low-level function. It does not perform any verification related to the
    inputs.

# Returns

- `T`: Field derivative with respect to `r`: `∂V/∂r`.
- `T`: Field derivative with respect to `ϕ`: `∂V/∂ϕ`.
- `T`: Field derivative with respect to `θ`: `∂V/∂θ`.
"""
function _igrf_geomagnetic_potential_gradient(
    n_max::Int,
    idx::Int,
    r_km::T,
    θ::T,
    ϕ::T,
    Δt::T,
    extrapolate::Bool,
    P::AbstractMatrix,
    dP::AbstractMatrix,
) where T<:Number
    # Auxiliary variables to select the IGRF coefficients.
    a = T(_IGRF_A)
    H = _IGRF_H
    G = _IGRF_G

    # Auxiliary variables to improve computational speed.
    sin_ϕ,  cos_ϕ  = sincos(ϕ)
    ratio = a / r_km
    fact  = ratio

    # Initialization of variables
    # ======================================================================================

    dVr = T(0)  # ........................ Derivative of the Geomagnetic potential w.r.t. r.
    dVθ = T(0)  # ........................ Derivative of the Geomagnetic potential w.r.t. θ.
    dVϕ = T(0)  # ........................ Derivative of the Geomagnetic potential w.r.t. ϕ.
    ΔG  = T(0)  # .................... Auxiliary variable to interpolate the G coefficients.
    ΔH  = T(0)  # .................... Auxiliary variable to interpolate the H coefficients.
    kg  = 1     # ............................ Index to obtain the values of the matrix `G`.
    kh  = 1     # ............................ Index to obtain the values of the matrix `H`.

    # Compute the Schmidt quasi-normalized associated Legendre functions and their first
    # order derivative, neglecting the phase term.
    legendre!(Val(:schmidt), P, θ, n_max, n_max; ph_term = false)
    dlegendre!(Val(:schmidt), dP, θ, P, n_max, n_max; ph_term = false)

    @inbounds for n in 1:n_max
        aux_dVr = T(0)
        aux_dVθ = T(0)
        aux_dVϕ = T(0)

        # Compute the contributions when `m = 0`
        # ==================================================================================

        # Get the coefficients in the epoch and interpolate to the desired time.
        Gnm_e0 = T(G[kg, idx])

        # If we need to extrapolate beyond the date we have measurements, we will use the
        # last column of the coefficients matrices that contains the expected
        # time-derivative.
        if !extrapolate
            Gnm_e1 = T(G[kg, idx+1])
            ΔG     = (Gnm_e1 - Gnm_e0) / T(5)
        else
            ΔG = T(G[kg,end])
        end

        Gnm  = Gnm_e0 + ΔG * Δt
        kg  += 1

        aux_dVr += -(n + 1) / r_km * Gnm * P[n+1, 1]
        aux_dVθ += Gnm * dP[n+1, 1]

        # Sine and cosine with m = 1
        # ==================================================================================
        #
        # This values will be used to update recursively `sin(m * ϕ)` and `cos(m * ϕ)`,
        # reducing the computational burden.
        #
        # TODO: Cache the computation.
        # We tried to compute those values only once using an external vector to store the
        # values. However, it leads to a worst performance. This behavior need further
        # investigation.
        sin_mϕ   = +sin_ϕ    # sin( 1 * λ_gc)
        sin_m_1ϕ = T(0)      # sin( 0 * λ_gc)
        sin_m_2ϕ = -sin_ϕ    # sin(-1 * λ_gc)
        cos_mϕ   = +cos_ϕ    # cos( 1 * λ_gc)
        cos_m_1ϕ = T(1)      # cos( 0 * λ_gc)
        cos_m_2ϕ = +cos_ϕ    # cos(-2 * λ_gc)

        # Other auxiliary variables that depend only on `n`
        # ==================================================================================

        fact_dVr = T(n + 1) / r_km

        # Compute the contributions when `m ∈ [1, n]`
        # ==================================================================================

        for m in 1:n
            # Compute recursively `sin(m * ϕ)` and `cos(m * ϕ)`.
            sin_mϕ = 2cos_ϕ * sin_m_1ϕ - sin_m_2ϕ
            cos_mϕ = 2cos_ϕ * cos_m_1ϕ - cos_m_2ϕ

            # Compute the coefficients `G_nm` and `H_nm`
            # ==============================================================================

            # Get the coefficients in the epoch and interpolate to the desired time.
            Gnm_e0 = T(G[kg, idx])
            Hnm_e0 = T(H[kh, idx])

            # If we need to extrapolate beyond the date we have measurements, we will use
            # the last column of the coefficients matrices that contains the expected
            # time-derivative.
            if !extrapolate
                Gnm_e1 = T(G[kg, idx+1])
                Hnm_e1 = T(H[kh, idx+1])
                ΔG     = (Gnm_e1 - Gnm_e0) / T(5)
                ΔH     = (Hnm_e1 - Hnm_e0) / T(5)
            else
                ΔG = T(G[kg, end])
                ΔH = T(H[kh, end])
            end

            Gnm  = Gnm_e0 + ΔG * Δt
            Hnm  = Hnm_e0 + ΔH * Δt
            kg  += 1
            kh  += 1

            GcHs_nm = Gnm * cos_mϕ + Hnm * sin_mϕ
            GsHc_nm = Gnm * sin_mϕ - Hnm * cos_mϕ

            # Compute the contributions for `m`
            # ==============================================================================

            aux_dVr += -fact_dVr * GcHs_nm * P[n+1, m+1]
            aux_dVθ += GcHs_nm * dP[n+1, m+1]
            aux_dVϕ += (θ == 0) ? -m * GsHc_nm * dP[n+1, m+1] : -m * GsHc_nm * P[n+1, m+1]

            # Update the values for the next step
            # ==============================================================================

            sin_m_2ϕ = sin_m_1ϕ
            sin_m_1ϕ = sin_mϕ
            cos_m_2ϕ = cos_m_1ϕ
            cos_m_1ϕ = cos_mϕ
        end

        # Perform final computations related to the summation in `n`
        # ==================================================================================

        # fact = (a / r)^(n + 1)
        fact *= ratio

        # aux_<> *= (a / r)^(n + 1)
        aux_dVr *= fact
        aux_dVϕ *= fact
        aux_dVθ *= fact

        dVr += aux_dVr
        dVϕ += aux_dVϕ
        dVθ += aux_dVθ
    end

    dVr *= a
    dVϕ *= a
    dVθ *= a

    return dVr, dVϕ, dVθ
end
