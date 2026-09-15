# AGENTS.md

Guidance for coding agents working on SatelliteToolboxGeomagneticField.jl, a Julia package providing geomagnetic field models (IGRF v14 and a simplified dipole model). Part of the [SatelliteToolbox.jl](https://github.com/JuliaSpace/SatelliteToolbox.jl) ecosystem.

## Package Structure

- Supported Julia versions: 1.10 through 1.13 (`[compat] julia = "1.10, 1.11, 1.12, 1.13"` — this list excludes future minors until updated).
- Entrypoint `src/SatelliteToolboxGeomagneticField.jl` includes, in order: `dipole/dipole_coefficients.jl`, `igrf/igrf_coefficients.jl` (constants), then `dipole/dipole.jl`, `igrf/igrf.jl` (functions), and ends with a PrecompileTools workload. New code can only reference symbols defined in earlier includes.
- `ext/SatelliteToolboxGeomagneticFieldZygoteExt.jl` is a package extension that defines a `ChainRulesCore.rrule` for `igrf` using ForwardDiff. It loads only when **both** ForwardDiff and Zygote are loaded (`[weakdeps]`/`[extensions]` in `Project.toml`), so exercising or testing it requires both packages.
- Test-only dependencies are declared via `[extras]` + `[targets]` in `Project.toml` (Test, DelimitedFiles, Aqua, AllocCheck, JET, ForwardDiff, Zygote), each with a `[compat]` entry so that Aqua's compat check passes. There is no `test/Project.toml`.

## Commands

- Instantiate: `julia --project=. -e 'using Pkg; Pkg.instantiate()'`
- Full test suite: `julia --project=. -e 'using Pkg; Pkg.test()'`
- Focused test file (run from the repo root; the `cd("test")` is required because the IGRF tests load data files via paths relative to the working directory): `julia --project=. -e 'using Test, DelimitedFiles, LinearAlgebra, ReferenceFrameRotations, StaticArrays, SatelliteToolboxGeomagneticField; using SatelliteToolboxBase: LowerTriangularStorage, RowMajor; cd("test"); include("igrf.jl")'` (same recipe with `include("dipole.jl")` for the dipole tests). Caveat: DelimitedFiles is a test-only dependency, so this works only if it is installed in your default environment (found via environment stacking); otherwise install it there first.
- Format: `julia -e 'using JuliaFormatter; format(".")'` — no `--project=.`; JuliaFormatter is not a package dependency and must be available in your default environment. Check for unwanted diffs with `git diff` afterward.
- Use generous timeouts: the first `Pkg.instantiate()`/`Pkg.test()` triggers precompilation and can run for minutes while printing little. Slow startup is precompilation, not a hang.
- There is no test-name selector; run a whole test file or the full suite.

## Tests

- `test/runtests.jl` always runs `test/igrf.jl` and `test/dipole.jl`, mirroring the two models in `src/`.
- On release Julia versions (not prereleases), `runtests.jl` then loads Aqua, AllocCheck, JET, ForwardDiff, and Zygote and runs `test/performance.jl` (Aqua, JET, allocation checks) and `test/zygote_extension.jl`. On prerelease/nightly Julia these are skipped (the packages are still installed by `Pkg.test`, so a resolver failure on nightly means one of them does not support that Julia yet).
- There are no platform skips: JET and AllocCheck run on every release Julia version, including 1.12+ on macOS.
- The IGRF reference values in `test/IGRF14_test_*.txt` were generated with the official `igrf14syn` Fortran routine (see https://github.com/JuliaSpace/IGRF_Test) and are compared with `atol = 3e-1` nT. Do not regenerate or loosen them to make tests pass.
- The allocation checks are a hard requirement: `igrf`/`igrfd` called with preallocated `P`/`dP` matrices and `show_warnings = Val(false)` must remain allocation-free (this is why the warning is controlled by a compile-time `Val` flag — `@warn` allocates, and AllocCheck flags it even in a branch that is not taken at runtime). The `MVector`s used inside the IGRF kernel are stack allocated only because they do not escape the function; do not return or store them.

## Code Style

- JuliaFormatter with `style = "blue"` plus the alignment options in `.JuliaFormatter.toml` (source of truth). CI does **not** run a format check, so format before committing.
- 92-character line limit. Comments are complete sentences ending with a period. Files start with a `## Description ###...` header block; sections are separated by `# == Section Name ==...` comment rules.
- Every function, macro, and struct has a docstring; keyword arguments are documented with their defaults as `(**Default** = ...)`.
- The coefficient matrices in `*_coefficients.jl` are wrapped in `#! format: off` — never reformat or hand-edit their values; they come verbatim from the official IGRF distribution (`igrf14coeffs.txt`) and the Kyoto WDC pole tables. Constants that depend on the tables (`_IGRF_NUM_EPOCHS`, `_GEOMAGNETIC_DIPOLE_MODEL_YEARS`) are derived from them in the same files instead of being hard-coded.
- Commit messages: gitmoji shortcode + imperative summary of 50 characters or fewer (e.g. `:bug: Fix ...`, `:books: Update ...`), body wrapped at 72 characters.
- `CHANGELOG.md` is maintained manually with badge-style entries (`![Bugfix][badge-bugfix] ...`); add entries when preparing a release.

## Behavioral Constraints

- Public API: `igrf`, `igrfd`, and `geomagnetic_dipole_field`. Output element types are documented as the float-promoted input types (e.g. Float32 inputs yield Float32 output) — preserve this when changing promotions. `date` is deliberately excluded from the promotion. The geodetic `igrf` method must convert the output of `geodetic_to_geocentric` back to the promoted type, since SatelliteToolboxTransformations promotes it with the ellipsoid parameter type.
- The warning keyword of `igrf`/`igrfd` is `show_warnings::Val{Bool}` (`Val(true)`/`Val(false)`), not a `Bool`; the Zygote extension must keep the same default as `igrf`.
- `igrfd` has a single method with the representation `R` as an optional positional argument, and both `igrf`/`igrfd` forward `kwargs...`; keep the keyword list documented in the docstrings since it is not visible in the signatures.
- `geomagnetic_dipole_field` throws `DimensionMismatch` for position vectors without exactly three elements (the input is converted to an `SVector{3}` first).
- `igrf` accepts dates in `[1900, 2035]`; results after 2030 are extrapolations and emit a warning. These bounds are tied to the IGRF v14 coefficient tables — do not change them independently of a coefficient update.
- The Zygote `rrule` packs all four `igrf` inputs into one vector before calling `ForwardDiff.jacobian`; this is intentional (differentiating w.r.t. `date` alone fails because the internal type promotion excludes it).

## Not Configured

- No documentation build: `docs/` contains only README assets; there is no `docs/make.jl` or Documenter setup.
- No `Manifest.toml` is committed, no `deps/build.jl` (so `Pkg.build()` is a no-op), no pre-commit hooks, and no CI formatting or docs jobs.
- CI (`.github/workflows/`): `ci.yml` tests Julia 1.10 (oldest supported) and the latest stable 1.x on Ubuntu (x64), macOS (arm64), and Windows (x64), building then testing, with coverage uploaded to Codecov; `ci-nightly.yml` runs the same matrix on Julia nightly without coverage; CompatHelper and TagBot handle dependency bumps and releases.
