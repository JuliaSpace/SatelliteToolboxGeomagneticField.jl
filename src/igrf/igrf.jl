## Description #############################################################################
#
# International Geomagnetic Field Model.
#
## References ##############################################################################
#
# [1] https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
# [2] https://www.ngdc.noaa.gov/IAGA/vmod/igrf12.f
# [3] https://www.mathworks.com/matlabcentral/fileexchange/34388-international-geomagnetic-reference-field--igrf--model
#
############################################################################################

export igrf, igrfd

############################################################################################
#                                        Functions                                         #
############################################################################################

"""
    igrfd(date::Number, <r, h>::T1, λ::T2, Ω::T3[, R]; kwargs...) where {T1<:Number, T2<:Number, T3<:Number} -> SVector{3, T}

**IGRF Model**

*Current version: v14*

Compute the geomagnetic field vector [nT] at the date `date` [Year A.D.] and position (`r`
or `h`, `λ`, `Ω`).

The position representation is defined by `R`. If `R` is `Val(:geocentric)`, the input must
be **geocentric** coordinates:

1. Distance from the Earth center `r` [m];
2. Geocentric latitude `λ` ∈ [-90°, +90°]; and
3. Geocentric longitude `Ω` ∈ [-180°, +180°].

If `R` is `Val(:geodetic)`, the input must be **geodetic** coordinates:

1. Altitude above the reference ellipsoid `h` (WGS-84) [m];
2. Geodetic latitude `λ` ∈ [-90°, +90°]; and
3. Geodetic longitude `Ω` ∈ [-180°, +180°].

If `R` is omitted, it defaults to `Val(:geocentric)`.

!!! warning

    We must have `1900 <= date <= 2035`. A warning message is printed for dates greater than
    2030 since the output is not reliable anymore. This message can be suppressed by setting
    the keyword `show_warnings` to `Val(false)`.

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
    matrices, it will be clamped. If it is equal to or lower than 0, it will be set to 1.
    (**Default** = 13)
- `show_warnings::Val`: If it is `Val(true)`, a warning is printed using `@warn` for dates
    greater than 2030, when the accuracy of the model is reduced. If it is `Val(false)`, the
    code related to the warning is removed at compile time, enabling allocation-free calls.
    (**Default** = `Val(true)`)
- `P::Union{Nothing, AbstractMatrix}`: An optional matrix that must contain at least
    `max_degree + 1 × max_degree + 1` real numbers that will be used to store the Legendre
    coefficients, reducing the allocations. If it is `nothing`, the matrix will be created
    when calling the function.
    (**Default** = `nothing`)
- `dP::Union{Nothing, AbstractMatrix}`: An optional matrix that must contain at least
    `max_degree + 1 × max_degree + 1` real numbers that will be used to store the Legendre
    derivative coefficients, reducing the allocations. If it is `nothing`, the matrix will
    be created when calling the function.
    (**Default** = `nothing`)

# Returns

- `SVector{3, T}`: Geomagnetic field vector [nT] at the desired location represented in the
    same input reference (geocentric or geodetic).

!!! info

    The output type `T` is obtained by promoting `T1`, `T2`, and `T3` to a float.
"""
function igrfd(
    date::Number,
    r::T1,
    λ::T2,
    Ω::T3,
    R::Union{Val{:geocentric}, Val{:geodetic}} = Val(:geocentric);
    kwargs...,
) where {T1 <: Number, T2 <: Number, T3 <: Number}
    T = promote_type(T1, T2, T3) |> float

    _check_latitude_and_longitude(λ, Ω, Val(:deg))

    return igrf(date, T(r), deg2rad(T(λ)), deg2rad(T(Ω)), R; kwargs...)
end

"""
    igrf(date::Number, <r, h>::T1, λ::T2, Ω::T3[, R]; kwargs...) where {T1<:Number, T2<:Number, T3<:Number} -> SVector{3, T}

**IGRF Model**

*Current version: v14*

Compute the geomagnetic field vector [nT] at the date `date` [Year A.D.] and position (`r`
or `h`, `λ`, `Ω`).

The position representation is defined by `R`. If `R` is `Val(:geocentric)`, the input must
be **geocentric** coordinates:

1. Distance from the Earth center `r` [m];
2. Geocentric latitude `λ` ∈ [-π / 2, +π / 2] [rad]; and
3. Geocentric longitude `Ω` ∈ [-π, +π] [rad].

If `R` is `Val(:geodetic)`, the input must be **geodetic** coordinates:

1. Altitude above the reference ellipsoid `h` (WGS-84) [m];
2. Geodetic latitude `λ` ∈ [-π / 2, +π / 2] [rad]; and
3. Geodetic longitude `Ω` ∈ [-π, +π] [rad].

If `R` is omitted, it defaults to `Val(:geocentric)`.

!!! warning

    We must have `1900 <= date <= 2035`. A warning message is printed for dates greater than
    2030 since the output is not reliable anymore. This message can be suppressed by setting
    the keyword `show_warnings` to `Val(false)`.

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
    matrices, it will be clamped. If it is equal to or lower than 0, it will be set to 1.
    (**Default** = 13)
- `show_warnings::Val`: If it is `Val(true)`, a warning is printed using `@warn` for dates
    greater than 2030, when the accuracy of the model is reduced. If it is `Val(false)`, the
    code related to the warning is removed at compile time, enabling allocation-free calls.
    (**Default** = `Val(true)`)
- `P::Union{Nothing, AbstractMatrix}`: An optional matrix that must contain at least
    `max_degree + 1 × max_degree + 1` real numbers that will be used to store the Legendre
    coefficients, reducing the allocations. If it is `nothing`, the matrix will be created
    when calling the function.
    (**Default** = `nothing`)
- `dP::Union{Nothing, AbstractMatrix}`: An optional matrix that must contain at least
    `max_degree + 1 × max_degree + 1` real numbers that will be used to store the Legendre
    derivative coefficients, reducing the allocations. If it is `nothing`, the matrix will
    be created when calling the function.
    (**Default** = `nothing`)

# Returns

- `SVector{3, T}`: Geomagnetic field vector [nT] at the desired location represented in the
    same input reference (geocentric or geodetic).

!!! info

    The output type `T` is obtained by promoting `T1`, `T2`, and `T3` to a float.
"""
function igrf(date::Number, r::Number, λ::Number, Ω::Number; kwargs...)
    return igrf(date, r, λ, Ω, Val(:geocentric); kwargs...)
end

function igrf(
    date::Number,
    r::T1,
    λ::T2,
    Ω::T3,
    ::Val{:geocentric};
    max_degree::Int = _IGRF_MAX_DEGREE,
    show_warnings::Val{S} = Val(true),
    P::Union{Nothing, AbstractMatrix} = nothing,
    dP::Union{Nothing, AbstractMatrix} = nothing,
) where {T1 <: Number, T2 <: Number, T3 <: Number, S}
    T = promote_type(T1, T2, T3) |> float

    # == Input Verification ================================================================

    # Check the date, since this model is valid for years between 1900 and
    # `_IGRF_LAST_YEAR`.
    if (date < 1900) || (date > _IGRF_LAST_YEAR)
        throw(
            ArgumentError(
                "This IGRF version will not work for years outside the interval [1900, $_IGRF_LAST_YEAR].",
            ),
        )
    end

    _check_latitude_and_longitude(λ, Ω, Val(:rad))

    # Warn the user that for dates after the year `_IGRF_RELIABLE_YEAR` the accuracy may be
    # reduced.
    # Since `S` is a compile-time constant, the code inside this branch is removed when
    # `show_warnings` is `Val(false)`, allowing allocation-free calls.
    if S && (date > _IGRF_RELIABLE_YEAR)
        @warn "The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than $_IGRF_RELIABLE_YEAR."
    end

    # == Input Variables Conversion ========================================================

    # Convert latitude / longitude to co-latitude and east-longitude.
    θ = T(π / 2) - T(λ)
    ϕ = T(Ω)

    # Check if the position is at one of the geographic poles, where the east component of
    # the field must be obtained by a limit since `sin(θ) = 0`.
    at_pole = (θ == 0) || (θ == T(π))

    # The input variable `r` is in [m], but all the algorithm requires it to be in [km].
    r_km = T(r) / 1000

    # == Preliminary Setup =================================================================

    # Compute the index of the epoch used to obtain the coefficients, which are provided
    # every `_IGRF_EPOCH_INTERVAL` years. Between two epochs, the coefficients are
    # interpolated. After the last epoch, they are extrapolated using the secular variation.
    idx   = clamp(floor(Int, (date - 1900) / _IGRF_EPOCH_INTERVAL) + 1, 1, _IGRF_NUM_EPOCHS)
    epoch = 1900 + (idx - 1) * _IGRF_EPOCH_INTERVAL

    # We must jump the first two columns that are reserved for the degree and order.
    idx += 2

    # Compute the elapsed time from the epoch of the coefficient selected by `idx`.
    Δt = T(date - epoch)

    # Compute the maximum spherical harmonic degree for the selected date, which is 10 for
    # the epochs before 1995, and clamp it with the one requested by the user.
    n_max = clamp(max_degree, 1, (epoch < 1995) ? 10 : _IGRF_MAX_DEGREE)

    # Check if the matrices related to Legendre must be computed.
    if isnothing(P)
        P = LowerTriangularStorage{RowMajor, T}(n_max + 1)

    else
        # If the user passed a matrix, we must check if there are enough space to store the
        # coefficients.
        rows, cols = size(P)

        if (rows < n_max + 1) || (cols < n_max + 1)
            throw(
                ArgumentError(
                    "Matrix `P` must have at least $(n_max + 1) rows and columns."
                ),
            )
        end
    end

    if isnothing(dP)
        dP = LowerTriangularStorage{RowMajor, T}(n_max + 1)

    else
        # If the user passed a matrix, we must check if there are enough space to store the
        # coefficients.
        rows, cols = size(dP)

        if (rows < n_max + 1) || (cols < n_max + 1)
            throw(
                ArgumentError(
                    "Matrix `dP` must have at least $(n_max + 1) rows and columns."
                ),
            )
        end
    end

    # == Geomagnetic Potential Gradient ====================================================

    dVr, dVϕ, dVθ = _igrf_geomagnetic_potential_gradient(
        n_max,
        idx,
        r_km,
        θ,
        ϕ,
        Δt,
        date >= _IGRF_LAST_YEAR_WITH_MEASUREMENTS,
        at_pole,
        P,
        dP,
    )

    # == Compute the Geomagnetic Field Vector in the Geocentric Reference Frame ============

    x = +dVθ / r_km
    y = at_pole ? -dVϕ / r_km : -dVϕ / (r_km * sin(θ))
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
    show_warnings::Val{S} = Val(true),
    P::Union{Nothing, AbstractMatrix} = nothing,
    dP::Union{Nothing, AbstractMatrix} = nothing,
) where {S}

    # TODO: This method has a small error (≈ 0.01 nT) compared with the `igrf12syn`.
    # However, the result is exactly the same as the MATLAB function in [3]. Hence, this
    # does not seem to be an error in the conversion from geodetic to geocentric
    # coordinates. This is probably caused by a numerical error. Further verification is
    # necessary.

    T = promote_type(typeof(h), typeof(λ), typeof(Ω)) |> float

    _check_latitude_and_longitude(λ, Ω, Val(:rad))

    # Convert the geodetic coordinates to geocentric coordinates. The conversion can promote
    # the values to the type of the ellipsoid parameters. Hence, we must convert the result
    # back to `T` to keep the documented output type.
    λ_gc, r = geodetic_to_geocentric(T(λ), T(h))

    # Compute the geomagnetic field in geocentric coordinates.
    B_gc = igrf(
        date,
        T(r),
        T(λ_gc),
        T(Ω),
        Val(:geocentric);
        max_degree = max_degree,
        show_warnings = show_warnings,
        P = P,
        dP = dP,
    )

    # Convert to geodetic coordinates.
    D_gd_gc = angle_to_dcm(T(λ_gc) - T(λ), :Y)
    B_gd    = D_gd_gc * B_gc

    return B_gd
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _check_latitude_and_longitude(λ::Number, Ω::Number, ::Val{:rad}) -> Nothing
    _check_latitude_and_longitude(λ::Number, Ω::Number, ::Val{:deg}) -> Nothing

Throw an `ArgumentError` if the latitude `λ` is outside the interval [-π / 2, +π / 2] rad or
if the longitude `Ω` is outside the interval [-π, +π] rad. If the last argument is
`Val(:deg)`, the same verification is performed considering that the inputs are in degrees.

The limits are converted to the floating-point type of the inputs before the comparison so
that, e.g., `Float32(π) / 2` is accepted as a valid Float32 latitude.
"""
function _check_latitude_and_longitude(λ::Number, Ω::Number, ::Val{:rad})
    if abs(λ) > oftype(float(λ), π) / 2
        throw(ArgumentError("The latitude must be between -π / 2 and +π / 2 rad."))
    end

    if abs(Ω) > oftype(float(Ω), π)
        throw(ArgumentError("The longitude must be between -π and +π rad."))
    end

    return nothing
end

function _check_latitude_and_longitude(λ::Number, Ω::Number, ::Val{:deg})
    if abs(λ) > 90
        throw(ArgumentError("The latitude must be between -90° and +90°."))
    end

    if abs(Ω) > 180
        throw(ArgumentError("The longitude must be between -180° and +180°."))
    end

    return nothing
end

"""
    _igrf_geomagnetic_potential_gradient(n_max::Int, idx::Int, r_km::T, θ::T, ϕ::T, Δt::T, extrapolate::Bool, at_pole::Bool, P::AbstractMatrix, dP::AbstractMatrix) where T<:Number -> NTuple{3, T}

Compute the gradient of the geomagnetic potential.

# Arguments

- `n_max::Int`: Maximum degree when computing the potential in the spherical harmonics.
- `idx::Int`: Column of the matrices `_IGRF_G` and `_IGRF_H` related to the desired epoch.
- `r_km::T`: Distance from the Earth's center [km].
- `θ::T`: Geocentric co-latitude [rad] ∈ [0, π].
- `ϕ::T`: East-longitude [rad].
- `Δt::T`: Elapsed time from the epoch related to the column `idx` [year].
- `extrapolate::Bool`: If `true`, the desired date is after the last epoch with
    measurements (`_IGRF_LAST_YEAR_WITH_MEASUREMENTS`). Hence, the coefficients are
    extrapolated using their time-derivative stored in the last column of the matrices
    `_IGRF_G` and `_IGRF_H`. Otherwise, they are interpolated between the columns `idx`
    and `idx + 1`.
- `at_pole::Bool`: If `true`, the position is at one of the geographic poles (`θ = 0` or
    `θ = π`). In this case, the derivative with respect to `ϕ` is replaced by the limit
    `∂V/∂ϕ / sin(θ)` so that the east component of the field can be obtained without
    dividing by `sin(θ) = 0`.
- `P::AbstractMatrix`: An auxiliary matrix to compute the values of the Legendre associated
    functions. It must have a dimension equal to or greater than `n_max + 1 × n_max + 1`.
- `dP::AbstractMatrix`: An auxiliary matrix to compute the derivatives of the Legendre
    associated functions. It must have a dimension equal to or greater than
    `n_max + 1 × n_max + 1`.

!!! warning

    This is a low-level function. It does not perform any verification related to the
    inputs.

# Returns

- `T`: Potential derivative with respect to `r`: `∂V/∂r`.
- `T`: Potential derivative with respect to `ϕ`: `∂V/∂ϕ`. If `at_pole` is `true`, the
    returned value is the limit of `∂V/∂ϕ / sin(θ)` at the pole instead.
- `T`: Potential derivative with respect to `θ`: `∂V/∂θ`.
"""
function _igrf_geomagnetic_potential_gradient(
    n_max::Int,
    idx::Int,
    r_km::T,
    θ::T,
    ϕ::T,
    Δt::T,
    extrapolate::Bool,
    at_pole::Bool,
    P::AbstractMatrix,
    dP::AbstractMatrix,
) where {T <: Number}
    # Auxiliary variables to select the IGRF coefficients.
    a = T(_IGRF_A)
    H = _IGRF_H
    G = _IGRF_G

    # Auxiliary variables to improve computational speed.
    sin_ϕ, cos_ϕ = sincos(ϕ)
    cos_θ = cos(θ)
    ratio = a / r_km
    fact = ratio

    # == Linear Model of the Coefficients ==================================================
    #
    # The coefficients at the desired date are obtained by the linear combination:
    #
    #   Cnm = w₀ * C[k, c₀] + w₁ * C[k, c₁],
    #
    # where `C` is one of the coefficient matrices. If the date lies between two epochs,
    # the columns are those epochs and the weights interpolate them. If the date is after
    # the last epoch with measurements, the second column is the secular variation and the
    # weights extrapolate the coefficients linearly, as in the reference implementation.
    # Both matrices have the same number of columns.
    c₀ = idx

    if extrapolate
        c₁ = size(G, 2)
        w₀ = T(1)
        w₁ = Δt
    else
        c₁ = idx + 1
        w₁ = Δt / T(_IGRF_EPOCH_INTERVAL)
        w₀ = T(1) - w₁
    end

    # == Sine and Cosine of the Multiples of the Longitude =================================
    #
    # Compute `sin(m * ϕ)` and `cos(m * ϕ)` for `m ∈ [0, n_max]` once using the Chebyshev
    # recurrence. The vectors are stack allocated since they do not escape this function.
    sin_mϕ = MVector{_IGRF_MAX_DEGREE + 1, T}(undef)
    cos_mϕ = MVector{_IGRF_MAX_DEGREE + 1, T}(undef)

    @inbounds begin
        sin_mϕ[1] = T(0)
        cos_mϕ[1] = T(1)
        sin_mϕ[2] = sin_ϕ
        cos_mϕ[2] = cos_ϕ

        for m in 2:n_max
            sin_mϕ[m + 1] = 2cos_ϕ * sin_mϕ[m] - sin_mϕ[m - 1]
            cos_mϕ[m + 1] = 2cos_ϕ * cos_mϕ[m] - cos_mϕ[m - 1]
        end
    end

    # == Initialization of Variables =======================================================

    dVr = T(0)  # ........................ Derivative of the Geomagnetic potential w.r.t. r.
    dVθ = T(0)  # ........................ Derivative of the Geomagnetic potential w.r.t. θ.
    dVϕ = T(0)  # ........................ Derivative of the Geomagnetic potential w.r.t. ϕ.
    kg  = 1     # ............................ Index to obtain the values of the matrix `G`.
    kh  = 1     # ............................ Index to obtain the values of the matrix `H`.

    # Compute the Schmidt quasi-normalized associated Legendre functions and their first
    # order derivative, neglecting the phase term. We must pass the maximum degree and
    # order explicitly. Otherwise, those functions will infer them from the matrix
    # dimensions, leading to unnecessary computations or even errors if the user passes
    # matrices larger than `n_max + 1 × n_max + 1`.
    legendre!(P, θ, _IGRF_LEGENDRE_COEFFICIENTS, n_max, n_max; ph_term = false)
    dlegendre!(dP, θ, P, _IGRF_LEGENDRE_COEFFICIENTS, n_max, n_max; ph_term = false)

    @inbounds for n in 1:n_max
        aux_dVr = T(0)
        aux_dVθ = T(0)
        aux_dVϕ = T(0)

        # == Compute the Contributions When `m = 0` ========================================

        Gnm = w₀ * T(G[kg, c₀]) + w₁ * T(G[kg, c₁])
        kg += 1

        aux_dVr += Gnm * P[n + 1, 1]
        aux_dVθ += Gnm * dP[n + 1, 1]

        # == Compute the Contributions When `m ∈ [1, n]` ===================================

        for m in 1:n
            s_mϕ = sin_mϕ[m + 1]
            c_mϕ = cos_mϕ[m + 1]

            # == Compute the Coefficients `G_nm` and `H_nm` ================================

            Gnm = w₀ * T(G[kg, c₀]) + w₁ * T(G[kg, c₁])
            Hnm = w₀ * T(H[kh, c₀]) + w₁ * T(H[kh, c₁])
            kg += 1
            kh += 1

            GcHs_nm = Gnm * c_mϕ + Hnm * s_mϕ
            GsHc_nm = Gnm * s_mϕ - Hnm * c_mϕ

            # == Compute the Contributions for `m` =========================================

            P_nm  = P[n + 1, m + 1]
            dP_nm = dP[n + 1, m + 1]

            aux_dVr += GcHs_nm * P_nm
            aux_dVθ += GcHs_nm * dP_nm

            # At the poles, `P_nm / sin(θ)` tends to `dP_nm * cos(θ)` for `m = 1` and to 0
            # for `m > 1`, when `dP_nm` is also 0. Hence, we can use the derivative to obtain
            # the limit of the east component, as in the reference implementation.
            aux_dVϕ += at_pole ? -m * GsHc_nm * dP_nm * cos_θ : -m * GsHc_nm * P_nm
        end

        # == Perform Final Computations Related to the Summation in `n` ====================

        # fact = (a / r)^(n + 1)
        fact *= ratio

        # The derivative with respect to `r` of `(a / r)^(n + 1)` adds the factor
        # `-(n + 1) / r`. The division by `r` is performed only once after the loop.
        dVr += -(n + 1) * fact * aux_dVr
        dVϕ += fact * aux_dVϕ
        dVθ += fact * aux_dVθ
    end

    dVr *= a / r_km
    dVϕ *= a
    dVθ *= a

    return dVr, dVϕ, dVθ
end
