"""Helper functions to compute gradients and isothermas for the
`isa_atmosphere` script"""
from math import exp

from constants import GRAVITY, AIR_GAS_CONSTANT


def grad(pressure_1, density_1, temperature_1,
         gradient_slope, altitude_1, altitude):
    """Computes air properties for an altitude that lies within a gradient
    section with respect to the International Standard Atmosphere

    Parameters
    ----------
    pressure_1 : float
        pressure at the starting point of the checkpoint
    density_1 : float
        pressure at the starting point of the checkpoint
    temperature_1 : float
        pressure at the starting point of the checkpoint
    gradient_slope : float
        slope of the gradient that is being calculated
    altitude_1 : float
        pressure at the starting point of the checkpoint
    altitude : float
        altitude where the properties will be computed
    Returns
    -------
    temperature : float
        temperature at the altitude
    pressure : float
        pressure at the altitude
    density : float
        density at the altitude
    """
    temperature = temperature_1 + gradient_slope * (altitude - altitude_1)
    pressure = pressure_1 * (temperature / temperature_1) ** (
        -GRAVITY / (gradient_slope * AIR_GAS_CONSTANT)
    )
    density = density_1 * (temperature / temperature_1) ** -(
        (GRAVITY / (gradient_slope * AIR_GAS_CONSTANT)) + 1
    )
    return (temperature, pressure, density)


def iso(pressure_1, temperature, altitude, altitude_1, density_1):
    """Computes air properties for an altitude that lies within a isothermal
    section with respect to the International Standard Atmosphere

    Parameters
    ----------
    pressure_1 : float
        pressure at the starting point of the checkpoint
    temperature : float
        temperature of the isothermal
    altitude : float
        altitude where properties will be computed
    altitude_1 : float
        altitude at the starting point of the checkpoint
    density_1 : float
        density at the starting point of the checkpoint

    Returns
    -------
    pressure : float
        pressure at the altitude
    density : float
        density at the altitude
    """
    pressure = pressure_1 * exp(
        -(GRAVITY / (AIR_GAS_CONSTANT * temperature)) * (altitude - altitude_1)
    )
    density = density_1 * exp(
        -(GRAVITY / (AIR_GAS_CONSTANT * temperature)) * (altitude - altitude_1)
    )
    return (pressure, density)
