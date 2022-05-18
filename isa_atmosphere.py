"""
This script computes the properties of a given altitude in the ISA standard
regime, such properties include air's Temperature, Pressure and Density.
"""
from utils import grad, iso
from constants import (
    T1,
    T2,
    T3,
    A1,
    A2,
    A3,
    A4,
    H_ISO11,
    H_ISO12,
    H_ISO21,
    H_ISO22,
    H_ISO31,
    H_ISO32,
    H_TOP,
    TEMPERATURE_SEA_LEVEL,
    PRESSURE_SEA_LEVEL,
    DENSITY_SEA_LEVEL,
)

# Bottom of 1 Isotherma Data
_, p, rho = grad(
    PRESSURE_SEA_LEVEL, DENSITY_SEA_LEVEL,
    TEMPERATURE_SEA_LEVEL, A1, 0, H_ISO11
)
p_iso11 = p
rho_iso11 = rho

# Top of 1 Isotherma Data
p, rho = iso(p, T1, H_ISO12, H_ISO11, rho)
p_iso12 = p
rho_iso12 = rho

# Bottom of 2 Isotherma Data
_, p, rho = grad(p, rho, T1, A2, H_ISO12, H_ISO21)
p_iso21 = p
rho_iso21 = rho

# Top of 2 Isotherma Data
p, rho = iso(p, T2, H_ISO22, H_ISO21, rho)
p_iso22 = p
rho_iso22 = rho

# Bottom of 3 Isotherma Data
_, p, rho = grad(p, rho, T2, A3, H_ISO22, H_ISO31)
p_iso31 = p
rho_iso31 = rho

# Top of 3 isotherma Data
p, rho = iso(p, T3, H_ISO32, H_ISO31, rho)
p_iso32 = p
rho_iso32 = rho

# Top of ISA table Data
_, p, rho = grad(p, rho, T3, A4, H_ISO32, H_TOP)
p_top = p
rho_top = rho


def atmosphere(altitude: float):
    """Computes air properties based on the International Standard
    Atmosphere for a given altitude

    Parameters
    ----------
    altitude : float
        altitude above sea level in meters

    Returns
    -------
    temperature : float
        air temperature at the give altitude in Kelvin
    pressure : float
        air pressure at the given altitude in Pascals
    density : float
        air density at the given altitude in kg/m³
    """
    if altitude == 0:
        pressure = PRESSURE_SEA_LEVEL
        density = DENSITY_SEA_LEVEL
        temperature = TEMPERATURE_SEA_LEVEL
    elif altitude > H_ISO32:
        temperature, pressure, density = grad(
            p_iso32, rho_iso32, T3, A4, H_ISO32, altitude
        )
    elif H_ISO32 >= altitude >= H_ISO31:
        pressure, density = iso(p_iso31, T3, altitude, H_ISO31, rho_iso31)
        temperature = T3
    elif H_ISO31 > altitude > H_ISO22:
        temperature, pressure, density = grad(
            p_iso22, rho_iso22, T2, A3, H_ISO22, altitude
        )
    elif H_ISO22 >= altitude >= H_ISO21:
        pressure, density = iso(p_iso21, T2, altitude, H_ISO21, rho_iso21)
        temperature = T2
    elif H_ISO21 > altitude > H_ISO12:
        temperature, pressure, density = grad(
            p_iso12, rho_iso12, T1, A2, H_ISO12, altitude
        )
    elif H_ISO12 >= altitude >= H_ISO11:
        pressure, density = iso(p_iso11, T1, altitude, H_ISO11, rho_iso11)
        temperature = T1
    elif H_ISO11 > altitude > 0:
        temperature, pressure, density = grad(
            PRESSURE_SEA_LEVEL,
            DENSITY_SEA_LEVEL,
            TEMPERATURE_SEA_LEVEL,
            A1,
            0,
            altitude,
        )
    return temperature, pressure, density
