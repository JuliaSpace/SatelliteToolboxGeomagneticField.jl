SatelliteToolboxGeomagneticField.jl Changelog
=============================================

Version 1.1.1
-------------

- ![Bugfix][badge-bugfix] The coefficients for the dipole model were not updated given the
  new IGRF v14 model.

Verison 1.1.0
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

[badge-breaking]: https://img.shields.io/badge/BREAKING-red.svg
[badge-deprecation]: https://img.shields.io/badge/Deprecation-orange.svg
[badge-feature]: https://img.shields.io/badge/Feature-green.svg
[badge-enhancement]: https://img.shields.io/badge/Enhancement-blue.svg
[badge-bugfix]: https://img.shields.io/badge/Bugfix-purple.svg
[badge-info]: https://img.shields.io/badge/Info-gray.svg

[gh-issue-3]: https://github.com/JuliaSpace/SatelliteToolboxGeomagneticField.jl/issues/3
