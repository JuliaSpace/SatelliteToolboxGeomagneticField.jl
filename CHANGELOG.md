SatelliteToolboxGeomagneticField.jl Changelog
=============================================

Version 2.0.0
-------------

- ![BREAKING][badge-breaking] The keywords `show_warnings::Bool` and `verbose::Val` of
  `igrf` and `igrfd` were merged into `show_warnings::Val{Bool}` (**Default** =
  `Val(true)`). Use `show_warnings = Val(false)` to suppress the warning about the reduced
  accuracy for dates after 2030, which also removes the related code at compile time,
  keeping the calls allocation-free. The Zygote extension now uses the same default as
  `igrf`, whereas it previously silenced the warning.
- ![BREAKING][badge-breaking] `geomagnetic_dipole_field` throws a `DimensionMismatch` if the
  position vector does not have three elements. Before, a longer vector silently returned a
  wrong field.
- ![Enhancement][badge-enhancement] The methods of `igrfd` were collapsed into a single one
  that forwards the keywords to `igrf`, and the range checks of the latitude and longitude
  were unified. The public signatures did not change.
- ![Enhancement][badge-enhancement] The IGRF kernel obtains the coefficients as a linear
  combination of two columns of the coefficient matrices, as the reference implementation
  `igrf14syn` does, removing the interpolation and extrapolation branches from the inner
  loop. The number of epochs is now derived from the coefficient matrix.
- ![Enhancement][badge-enhancement] `igrf` and `igrfd` are about 16% faster with
  preallocated matrices since the sines and cosines of the multiples of the longitude are
  computed only once per call.
- ![Enhancement][badge-enhancement] `geomagnetic_dipole_field` is about 30% faster since the
  interval search in the pole table uses `searchsortedlast`.
- ![Bugfix][badge-bugfix] The geodetic methods of `igrf` and `igrfd` returned a `Float64`
  vector for `Float32` inputs when used with SatelliteToolboxTransformations.jl v1.3,
  contradicting the documented output type.
- ![Bugfix][badge-bugfix] The special case used to compute the east component of the field
  at the geographic poles only covered the north pole. The south pole is now handled with
  the same limit used by `igrf14syn`.
- ![Info][badge-info] The package now supports Julia 1.13 and requires
  SatelliteToolboxBase.jl v2 and SatelliteToolboxTransformations.jl v1.3.
- ![Info][badge-info] The test dependencies are declared in `Project.toml` instead of being
  installed when the test suite runs, and the JET and AllocCheck tests are no longer skipped
  on Julia 1.12+.
- ![Info][badge-info] The references now cite IGRF-14, the docstrings state that the output
  is represented in the NED reference system, and the README examples were regenerated.

Version 1.3.1
-------------

- ![Enhancement][badge-enhancement] The maximum degree and order are now passed explicitly
  to the associated Legendre functions. This modification fixes an error when the matrices
  `P` and `dP` were larger than `14 × 14` and improves the performance when `max_degree` is
  lower than 13.
- ![Enhancement][badge-enhancement] The package now has a precompilation workload, highly
  reducing the time to first call of `igrf`, `igrfd`, and `geomagnetic_dipole_field`.
- ![Bugfix][badge-bugfix] The function `geomagnetic_dipole_field` was not working with
  vectors of integers, although the documentation stated that the input element type is
  converted to a float.
- ![Bugfix][badge-bugfix] The geodetic methods of `igrf` and `igrfd` were not validating
  the input latitude and longitude, silently returning meaningless results for values
  outside the valid range.
- ![Bugfix][badge-bugfix] Some error messages contained wrong information about the input
  units and the accepted date interval.
- ![Info][badge-info] The documentation received several fixes, including outdated examples
  in README.md that called a function that no longer exists.

Version 1.3.0
-------------

- ![Feature][badge-feature] The package has now differentiability support. (PR
  [#4][gh-pr-4])

Version 1.2.0
-------------

- ![Enhancement][badge-enhancement] The IGRF algorithm uses by default the
  `LowerTriangularStorage` from SatelliteToolboxBase.jl for the matrices `P` and `dP` if
  they are not provided. This reduces the memory footprint by half without a noticeable
  impact on performance.

Version 1.1.1
-------------

- ![Bugfix][badge-bugfix] The coefficients for the dipole model were not updated given the
  new IGRF v14 model.

Version 1.1.0
-------------

- ![Feature][badge-feature] Update IGRF to v14. (Issue [#3][gh-issue-3])

Version 1.0.0
-------------

- ![Info][badge-info] We dropped support for Julia 1.6. This version only supports the
  current Julia version and v1.10 (LTS).
- ![Info][badge-info] This version does not have breaking changes. We bump the version to
  1.0.0 because we now consider the API stable.

Version 0.1.2
-------------

- ![Enhancement][badge-enhancement] Minor source-code updates.

Version 0.1.1
-------------

- ![Enhancement][badge-enhancement] Documentation update.

Version 0.1.0
-------------

- Initial version.
  - This version was based on the code in **SatelliteToolbox.jl**.

[badge-breaking]: https://img.shields.io/badge/Breaking-DC2626?style=flat-square
[badge-deprecation]: https://img.shields.io/badge/Deprecation-D97706?style=flat-square
[badge-feature]: https://img.shields.io/badge/Feature-16A34A?style=flat-square
[badge-enhancement]: https://img.shields.io/badge/Enhancement-0284C7?style=flat-square
[badge-bugfix]: https://img.shields.io/badge/Bugfix-DB2777?style=flat-square
[badge-info]: https://img.shields.io/badge/Info-475569?style=flat-square

[gh-issue-3]: https://github.com/JuliaSpace/SatelliteToolboxGeomagneticField.jl/issues/3

[gh-pr-4]: https://github.com/JuliaSpace/SatelliteToolboxGeomagneticField.jl/pull/4
