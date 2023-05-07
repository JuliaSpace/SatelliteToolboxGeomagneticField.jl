SatelliteToolboxGeomagneticField.jl
===================================

[![CI](https://github.com/JuliaSpace/SatelliteToolboxGeomagneticField.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaSpace/SatelliteToolboxGeomagneticField.jl/actions/workflows/ci.yml)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

This packages contains models to compute the geomagnetic field vector. We currently have two
models implemented:

1. The [International Geomagnetic Reference Field (IGRF) v13](https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html); and
2. The simplified dipole model.

## Installation

```julia
julia> using Pkg
julia> Pkg.install("SatelliteToolboxGeomagneticField")
```

## Usage

### IGRF

We have a native Julia implementation of the [International Geomagnetic Reference Field
(IGRF) v13](https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html) based on **[1]**. This mode can
be accessed by two functions: `irgf` and `irgfd`.

```julia
function igrfd(date::Number, <r, h>::T1, λ::T2, Ω::T3[, R]; kwargs...) where {T1<:Number, T2<:Number, T3<:Number}
```
Compute the geomagnetic field vector [nT] at the date `date` [Year A.D.] and position (`r` or `h`, `λ`, `Ω`).

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

> **Warning**
> We must have `1900 <= date <= 2030`. A warning message is printed for dates greater than
> 2025 since the output is not reliable anymore. This message can be suppressed by setting
> the keyword `show_warnings` to `false`.

> **Note**
> The output vector will be represented in the same reference system selected by the
> parameter `R` (geocentric or geodetic). The Y-axis of the output reference system always
> points East. In case of **geocentric coordinates**, the Z-axis points toward the center of
> Earth and the X-axis completes a right-handed coordinate system. In case of **geodetic
> coordinates**, the X-axis is tangent to the ellipsoid at the selected location and points
> toward North, whereas the Z-axis completes a right-hand coordinate system.

The following keywords are available:

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

This function returns a `SVector{3, T}`, which is the geomagnetic field vector [nT] at the
desired location represented in the same input reference (geocentric or geodetic). Notice
that the output type `T` is obtained by promoting `T1`, `T2`, and `T3` to a float.

``` julia
function igrf(date::Number, <r, h>::T1, λ::T2, Ω::T3[, R]; kwargs...) where {T1<:Number, T2<:Number, T3<:Number}
```

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

> **Warning**
> We must have `1900 <= date <= 2030`. A warning message is printed for dates greater than
> 2025 since the output is not reliable anymore. This message can be suppressed by setting
> the keyword `show_warnings` to `false`.

> **Note**
> The output vector will be represented in the same reference system selected by the
> parameter `R` (geocentric or geodetic). The Y-axis of the output reference system always
> points East. In case of **geocentric coordinates**, the Z-axis points toward the center of
> Earth and the X-axis completes a right-handed coordinate system. In case of **geodetic
> coordinates**, the X-axis is tangent to the ellipsoid at the selected location and points
> toward North, whereas the Z-axis completes a right-hand coordinate system.

The following keywords are available:

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

This function returns a `SVector{3, T}`, which is the geomagnetic field vector [nT] at the
desired location represented in the same input reference (geocentric or geodetic). Notice
that the output type `T` is obtained by promoting `T1`, `T2`, and `T3` to a float.

```julia
julia> igrf(2017.12313, 640e3, 50 * pi / 180, 25 * pi / 180, Val(:geodetic))
3-element SVector{3, Float64} with indices SOneTo(3):
 15365.787505205588
  1274.9958640696996
 34201.2182033379

julia> igrfd(2017.12313, 640e3, 50, 25, Val(:geodetic))
3-element SVector{3, Float64} with indices SOneTo(3):
 15365.787505205588
  1274.9958640696996
 34201.2182033379

julia> igrf(2017.12313, 6371e3 + 640e3, 50 * pi / 180, 25 * pi / 180, Val(:geocentric))
3-element SVector{3, Float64} with indices SOneTo(3):
 15165.486702524944
  1269.7264334427584
 34243.04928373083

julia> igrfd(2017.12313, 6371e3 + 640e3, 50, 25, Val(:geocentric))
3-element SVector{3, Float64} with indices SOneTo(3):
 15165.486702524944
  1269.7264334427584
 34243.04928373083
```

```julia
julia> igrf(2026, 6371e3 + 640e3, 50 * pi/180, 25 * pi/180)
┌ Warning: The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than 2025.
└ @ SatelliteToolboxGeomagneticField ~/tmp/SatelliteToolboxGeomagneticField/src/igrf/igrf.jl:301
3-element SVector{3, Float64} with indices SOneTo(3):
 15118.591511098812
  1588.129544718569
 34668.84185460438

julia> igrfd(2026, 6371e3 + 640e3, 50, 25)
┌ Warning: The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than 2025.
└ @ SatelliteToolboxGeomagneticField ~/tmp/SatelliteToolboxGeomagneticField/src/igrf/igrf.jl:301
3-element SVector{3, Float64} with indices SOneTo(3):
 15118.591511098812
  1588.129544718569
 34668.84185460438

julia> igrf(2026, 6371e3 + 640e3, 50 * pi / 180, 25 * pi / 180; show_warnings = false)
3-element SVector{3, Float64} with indices SOneTo(3):
 15118.591511098812
  1588.129544718569
 34668.84185460438

julia> igrfd(2026, 6371e3+640e3, 50, 25; show_warnings = false)
3-element SVector{3, Float64} with indices SOneTo(3):
 15118.591511098812
  1588.129544718569
 34668.84185460438
```

### Simplified dipole model

This model assumes that the Earth geomagnetic field is a perfect dipole. The approximation
is not good but it can be sufficient for some analysis, such as those carried out at the
Pre-Phase A of a space mission when the uncertainties are high.

The following function computes the geomagnetic field using the simplified dipole model:

```julia
function geomagnetic_dipole_field(r_e::AbstractVector{T}, year::Number = 2020) where T
```

Compute the geomagnetic field [nT] using the simplified dipole model at position `r_e` (ECEF
reference frame) [m]. This function uses the year `year` to obtain the position of the South
geomagnetic pole (which lies in the North hemisphere) and the dipole moment. If `year` is
omitted, it defaults to 2020.

Notice that:

1. The output vector will be represented in the ECEF reference frame.
2. The returned vector type is obtained by converting `T` to a float.
3. The south geomagnetic pole position and dipole moment is obtained by interpolating the
    values provided in **[2]**.
    
```julia
julia> r_e = [0; 0; R0 + 200e3];

julia> geomag_dipole(r_e)
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
   1286.0242861717802
  -4232.804339060699
 -53444.68086319672

julia> geomag_dipole(r_e, 1986)
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
   1715.2656071053527
  -4964.59806084178
 -54246.30480714959

julia> r_e = [R0+200e3;0;0];

julia> geomag_dipole(r_e)
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 -2572.0485723435604
 -4232.804339060699
 26722.34043159836

julia> geomag_dipole(r_e, 1986)
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 -3430.5312142107055
 -4964.59806084178
 27123.152403574793
```

## References

- **[1]** **Alken, P.; Thébault, E.; Beggan, C. D. et al. (2021)**. *International
  Geomagnetic Reference Field: The Thirteenth Generation*. **Earth Planets Space**, v. 73,
  n. 49.
- **[2]** [http://wdc.kugi.kyoto-u.ac.jp/poles/polesexp.html](http://wdc.kugi.kyoto-u.ac.jp/poles/polesexp.html)
