SatelliteToolboxGeomagneticField.jl Changelog
=============================================

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
