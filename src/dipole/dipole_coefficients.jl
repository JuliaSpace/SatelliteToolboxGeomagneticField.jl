## Description #############################################################################
#
# Coefficients for the dipole model for the Earth geomagnetic field.
#
## References ##############################################################################
#
# [1] http://wdc.kugi.kyoto-u.ac.jp/poles/polesexp.html
#
############################################################################################

# This matrix has the following format:
#   - 1st column: Year;
#   - 2nd column: South geomagnetic pole latitude [°];
#   - 3rd column: South geomagnetic pole longitude [°]; and
#   - 4th column: Dipole moment [10²² A.m²].
const _GEOMAGNETIC_DIPOLE_MODEL_COEFFICIENTS = [
    1900    +78.7   -68.8   +8.32
    1905    +78.7   -68.7   +8.30
    1910    +78.7   -68.7   +8.27
    1915    +78.6   -68.6   +8.24
    1920    +78.6   -68.4   +8.20
    1925    +78.6   -68.3   +8.16
    1930    +78.6   -68.3   +8.13
    1935    +78.6   -68.4   +8.11
    1940    +78.5   -68.5   +8.09
    1945    +78.5   -68.5   +8.08
    1950    +78.5   -68.8   +8.06
    1955    +78.5   -69.2   +8.05
    1960    +78.6   -69.5   +8.03
    1965    +78.6   -69.9   +8.00
    1970    +78.7   -70.2   +7.97
    1975    +78.8   -70.5   +7.94
    1980    +78.9   -70.8   +7.91
    1985    +79.0   -70.9   +7.87
    1990    +79.2   -71.1   +7.84
    1995    +79.4   -71.4   +7.81
    2000    +79.6   -71.6   +7.79
    2005    +79.8   -71.8   +7.77
    2010    +80.1   -72.2   +7.75
    2011    +80.1   -72.3   +7.74
    2012    +80.2   -72.4   +7.74
    2013    +80.3   -72.5   +7.73
    2014    +80.3   -72.5   +7.73
    2015    +80.4   -72.6   +7.72
    2016    +80.4   -72.6   +7.72
    2017    +80.5   -72.6   +7.72
    2018    +80.5   -72.7   +7.71
    2019    +80.6   -72.7   +7.71
    2020    +80.7   -72.7   +7.71
    2021    +80.7   -72.7   +7.71
    2022    +80.7   -72.7   +7.70
    2023    +80.8   -72.7   +7.70
    2024    +80.8   -72.7   +7.69
    2025    +80.8   -72.8   +7.69
    2026    +80.9   -72.8   +7.69
    2027    +80.9   -72.8   +7.68
    2028    +81.0   -72.9   +7.68
    2029    +81.0   -72.9   +7.67
    2030    +81.1   -73.0   +7.67
]
