# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Dipole model for the Earth geomagnetic field.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] http://helios.fmi.fi/~juusolal/geomagnetism/Lectures/Chapter3_dipole.pdf
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export geomagnetic_dipole_field

"""
    geomagnetic_dipole_field(r_e::AbstractVector{T}, year::Number = 2020) where T -> SVector{3, T}

Compute the geomagnetic field [nT] using the simplified dipole model at position `r_e` (ECEF
reference frame) [m]. This function uses the year `year` to obtain the position of the South
geomagnetic pole (which lies in the North hemisphere) and the dipole moment. If `year` is
omitted, it defaults to 2020.

# Remarks

1. The output vector will be represented in the ECEF reference frame.
2. The returned vector type is obtained by converting `T` to a float.
3. The south geomagnetic pole position and dipole moment is obtained by interpolating the
    values provided in **[1]**.

# References

- **[1]**: http://wdc.kugi.kyoto-u.ac.jp/poles/polesexp.html
"""
function geomagnetic_dipole_field(r_e::AbstractVector{T}, year::Number = 2020) where T
    # Obtain the geomagnetic dipole coefficients.
    pole_lat, pole_lon, m = _geomagnetic_dipole_coefficients(year::Number)

    # DCM that converts the ECEF into the geomagnetic coordinates.
    Dge = angle_to_dcm(T(pole_lon), T(π / 2) - T(pole_lat), :ZY)

    # Compute the dipole momentum represented in the ECEF reference frame.
    k₀_e = T(1e-7) * T(m) * (Dge' * SVector{3, T}(0, 0, -1))

    # Compute the distance from the Earth center of the desired point.
    r = norm(r_e)

    # Compute the unitary vector that points to the desired direction.
    er_e = SVector{3, T}(r_e) / r

    # Compute the geomagnetic field vector [nT].
    B_e = (3er_e * er_e' - I) * k₀_e * T(1e9) / r^3

    return B_e
end

############################################################################################
#                                    Private Functions
############################################################################################

# Obtain the geomagnetic dipole coefficients (latitude and longitude of the south
# geomagnetic pole and the dipole moment) in `year`.
function _geomagnetic_dipole_coefficients(year::Number)
    C = _GEOMAGNETIC_DIPOLE_MODEL_COEFFICIENTS

    # We will interpolate linearly the coefficients. Hence, we first check if the year is
    # beyond the limits, when we will clamp the output (flat extrapolation).
    if year <= C[1, 1]
        return deg2rad(C[1, 2]), deg2rad(C[1, 3]), C[1, 4] * 1e22

    elseif year >= C[end, 1]
        return deg2rad(C[end, 2]), deg2rad(C[end, 3]), C[end, 4] * 1e22

    else
        # Perform a interval binary search of `year` in `C[1, :]`. It means that this
        # algorithm returns `id` such that `C[id, 1] <= year < C[id + 1, 1]`.
        num_elements = size(C, 1)
        low  = 1
        high = num_elements
        id   = nothing

        while low < high
            mid = div(low + high, 2, RoundDown)

            if (mid == num_elements) || (C[mid, 1] <= year < C[mid + 1, 1])
                id = mid
                break

            elseif (C[mid, 1] < year)
                low = mid + 1

            elseif (C[mid, 1] > year)
                high = mid

            end
        end

        if isnothing(id)
            id = low
        end

        # Linearly interpolate the values.
        Δt    = year - C[id, 1]

        year₀ = C[id, 1]
        lat₀  = C[id, 2]
        lon₀  = C[id, 3]
        m₀    = C[id, 4]

        Δyear = C[id+1, 1] - year₀
        Δlat  = C[id+1, 2] - lat₀
        Δlon  = C[id+1, 3] - lon₀
        Δm    = C[id+1, 4] - m₀

        lat   = lat₀ + Δlat / Δyear * Δt
        lon   = lon₀ + Δlon / Δyear * Δt
        m     = m₀ + Δm / Δyear * Δt

        # Return the values with the correct units.
        return deg2rad(lat), deg2rad(lon), m * 1e22
    end
end
