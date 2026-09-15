## Description #############################################################################
#
# Dipole model for the Earth geomagnetic field.
#
## References ##############################################################################
#
# [1] http://helios.fmi.fi/~juusolal/geomagnetism/Lectures/Chapter3_dipole.pdf
#
############################################################################################

export geomagnetic_dipole_field

"""
    geomagnetic_dipole_field(r_e::AbstractVector{T}, year::Number = 2020) where T -> SVector{3, float(T)}

Compute the geomagnetic field [nT] using the simplified dipole model at position `r_e` (ECEF
reference frame) [m]. This function uses the year `year` to obtain the position of the South
geomagnetic pole (which lies in the North hemisphere) and the dipole moment. If `year` is
omitted, it defaults to 2020.

# Remarks

1. The output vector will be represented in the ECEF reference frame.
2. The returned vector type is obtained by converting `T` to a float.
3. The south geomagnetic pole position and dipole moment are obtained by interpolating the
    values provided in **[1]**.
4. A `DimensionMismatch` is thrown if `r_e` does not have three elements.

# References

- **[1]**: http://wdc.kugi.kyoto-u.ac.jp/poles/polesexp.html
"""
function geomagnetic_dipole_field(r_e::AbstractVector{T}, year::Number = 2020) where {T}
    # Convert the input type to a float to support, e.g., vectors of integers.
    Tf = float(T)

    # Obtain the geomagnetic dipole coefficients.
    pole_lat, pole_lon, m = _geomagnetic_dipole_coefficients(year)

    # DCM that converts the ECEF into the geomagnetic coordinates.
    Dge = angle_to_dcm(Tf(pole_lon), Tf(π / 2) - Tf(pole_lat), :ZY)

    # Compute the dipole momentum represented in the ECEF reference frame.
    k₀_e = Tf(1e-7) * Tf(m) * (Dge' * SVector{3, Tf}(0, 0, -1))

    # Convert the position to a static vector, which also throws a `DimensionMismatch` if
    # the input does not have three elements.
    r_e_s = SVector{3, Tf}(r_e)

    # Compute the distance from the Earth center of the desired point.
    r = norm(r_e_s)

    # Compute the unitary vector that points to the desired direction.
    er_e = r_e_s / r

    # Compute the geomagnetic field vector [nT].
    B_e = (3er_e * er_e' - I) * k₀_e * Tf(1e9) / r^3

    return B_e
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _geomagnetic_dipole_coefficients(year::Number) -> Float64, Float64, Float64

Obtain the geomagnetic dipole coefficients in `year` by linearly interpolating the table
`_GEOMAGNETIC_DIPOLE_MODEL_COEFFICIENTS`. If `year` is outside the interval of the table,
the values at the closest limit are returned (flat extrapolation).

# Returns

- `Float64`: Latitude of the south geomagnetic pole [rad].
- `Float64`: Longitude of the south geomagnetic pole [rad].
- `Float64`: Dipole moment [A.m²].
"""
function _geomagnetic_dipole_coefficients(year::Number)
    C     = _GEOMAGNETIC_DIPOLE_MODEL_COEFFICIENTS
    years = _GEOMAGNETIC_DIPOLE_MODEL_YEARS

    # Clamp the year to the limits of the table, leading to a flat extrapolation.
    if year <= years[begin]
        return deg2rad(C[begin, 2]), deg2rad(C[begin, 3]), C[begin, 4] * 1e22

    elseif year >= years[end]
        return deg2rad(C[end, 2]), deg2rad(C[end, 3]), C[end, 4] * 1e22
    end

    # Find `id` such that `years[id] <= year < years[id + 1]`, which exists since `year` is
    # strictly inside the interval of the table.
    id = searchsortedlast(years, year)

    # Linearly interpolate the values.
    year₀ = years[id]
    lat₀  = C[id, 2]
    lon₀  = C[id, 3]
    m₀    = C[id, 4]

    Δyear = years[id + 1] - year₀
    Δlat  = C[id + 1, 2] - lat₀
    Δlon  = C[id + 1, 3] - lon₀
    Δm    = C[id + 1, 4] - m₀

    Δt = (year - year₀) / Δyear

    lat = lat₀ + Δlat * Δt
    lon = lon₀ + Δlon * Δt
    m   = m₀ + Δm * Δt

    # Return the values with the correct units.
    return deg2rad(lat), deg2rad(lon), m * 1e22
end
