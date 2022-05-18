'''Script containing useful physical constants to be used for computations'''

GRAVITY = 9.81  # m/s2
AIR_GAS_CONSTANT = 287.05  # J/kgK

# Isothermal Temperatures
T1 = 216.66  # ºK, temperature at the top of the gradient 1
T2 = 282.66  # ºK, temperature at the top of the gradient 2
T3 = 165.66  # ºK, temperature at the top of the gradient 3

# Isotherma height intervals in meters (m)
H_ISO11 = 11e3
H_ISO12 = 25e3
H_ISO21 = 47e3
H_ISO22 = 53e3
H_ISO31 = 79e3
H_ISO32 = 90e3
H_TOP = 105e3  # altitude at the end of gradient 3

# Sea level properties
TEMPERATURE_SEA_LEVEL = 288.16  # ºK, temperature of air at sea level
PRESSURE_SEA_LEVEL = 1.01325e5  # N/m2, pressure of air at sea level
DENSITY_SEA_LEVEL = 1.225  # kg/m3, density of air at sea level

# Gradient values (slopes on ISA table)
A1 = -6.5e-3  # K/m, gradient 1
A2 = 3e-3  # K/m, gradient 2
A3 = -4.5e-3  # K/m, gradient 3
A4 = 4e-3  # K/m, gradient 4
