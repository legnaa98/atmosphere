# -*- coding: UTF8 -*-
from numpy import *
import matplotlib.pyplot as plt
'''
This script computes the properties of a given altitude in the ISA standard regime,
such properties include air's Temperature, Pressure and Density. The properties can be calculated for both
a single altitude or a set of altitudes within a given interval. This script also plots the curves for
T vs. altitude, P vs. altitude & Rho(density) vs. altitude.

Author: Alejandro Angel
Contact: alejandro.angelc@upb.edu.co
Date: january 23rd, 2019
'''


## Useful Constants
# ========================================================================================
# Isothermal Temperatures
T_t1 = 216.66 # ºK, temperature at the top of the gradient 1
T_t2 = 282.66 # ºK, temperature at the top of the gradient 2
T_t3 = 165.66 # ºK, temperature at the top of the gradient 3

# Isotherma height intervals in meters (m)
h_iso11 = 11e3
h_iso12 = 25e3
h_iso21 = 47e3
h_iso22 = 53e3
h_iso31 = 79e3
h_iso32 = 90e3
h_top = 105e3 # altitude at the end of gradient 3

# Physical Constants
g = 9.81 # m/s2, acceleration of gravity
R = 287.05 # J/kgK, air's gas constant

# Sea level properties
T_sl = 288.16 # ºK, temperature of air at sea level
p_sl = 1.01325e5 # N/m2, pressure of air at sea level
rho_sl = 1.225 # kg/m3, density of air at sea level

# Gradient values (slopes on ISA table)
a1 = -6.5e-3 # K/m, gradient 1
a2 = 3e-3 # K/m, gradient 2
a3 = -4.5e-3 # K/m, gradient 3
a4 = 4e-3 # K/m, gradient 4
# ========================================================================================
# Gradient and isotherma functions
def grad(p1,rho1,T1,a,h1,h):
	
	T = T1 + a*(h-h1)
	p = p1*(T/T1)**(-g/(a*R))
	rho = rho1*(T/T1)**-((g/(a*R))+1)
	return(T, p, rho)

def iso(p1,T,h,h1,rho1):

	p = p1*exp(-(g/(R*T))*(h-h1))
	rho = rho1*exp(-(g/(R*T))*(h-h1))
	return(p, rho)

# ========================================================================================
# Useful reference values

# Bottom of 1 Isotherma Data
T,p,rho = grad(p_sl,rho_sl,T_sl,a1,0,h_iso11)
p_iso11 = p; rho_iso11 = rho; #print(p_iso11,rho_iso11)

# Top of 1 Isotherma Data
p,rho = iso(p,T_t1,h_iso12,h_iso11,rho)
p_iso12 = p; rho_iso12 = rho; #print(p_iso12,rho_iso12)

# Bottom of 2 Isotherma Data
T,p,rho = grad(p,rho,T_t1,a2,h_iso12,h_iso21)
p_iso21 = p; rho_iso21 = rho; #print(p_iso21,rho_iso21)

# Top of 2 Isotherma Data
p,rho = iso(p,T_t2,h_iso22,h_iso21,rho)
p_iso22 = p; rho_iso22 = rho; #print(p_iso22,rho_iso22)

# Bottom of 3 Isotherma Data
T,p,rho = grad(p,rho,T_t2,a3,h_iso22,h_iso31)
p_iso31 = p; rho_iso31 = rho; #print(p_iso31,rho_iso31)

# Top of 3 isotherma Data
p,rho = iso(p,T_t3,h_iso32,h_iso31,rho)
p_iso32 = p; rho_iso32 = rho; #print(p_iso32,rho_iso32)

# Top of ISA table Data
T,p,rho = grad(p,rho,T_t3,a4,h_iso32,h_top)
p_top = p; rho_top = rho; #print(p_top,rho_top)
# ========================================================================================
usr_dir = input("Do you want to calculate the properties on a single altitude? (y/n): \n \t >> "); usr_dir = usr_dir.lower()

if usr_dir == "y":
	h = float(input("Enter the altitude in meters (m): ")); h = [h]

else:	
	h_b = float(input("Enter the lower altitude in meters (m): ")); h_t = float(input("Enter the higher altitude in meters (m): "))
	tot = int(input("Altitudes between the interval [" + str(round(h_b,2)) + ", "+ str(round(h_t,2)) + "] from which you'd like to compute their properties: "))
	h = linspace(h_b, h_t, tot)

print("\n")
# Empty lists to be filled with computed data
P = [] # pressure vector
RHO = [] # density vector
t = [] # temperature vector

# Test each altitude and compute its properties depending on its location on the ISA standard
for height in h:
	print("For a height of %f km \n" %(round(height/1000,2)))

	if height == 0: # Check wether the altitude corresponds to sea level
		p = p_sl
		rho = rho_sl
		T = T_sl
	
	# Test every interval from the top down

	elif height > h_iso32:
		T,p,rho = grad(p_iso32,rho_iso32,T_t3,a4,h_iso32,height)

	elif h_iso32 >= height >= h_iso31:
		p,rho = iso(p_iso31,T_t3,height,h_iso31,rho_iso31)
		T = T_t3

	elif h_iso31 > height > h_iso22:
		T,p,rho = grad(p_iso22,rho_iso22,T_t2,a3,h_iso22,height)

	elif h_iso22 >= height >= h_iso21:
		p,rho = iso(p_iso21,T_t2,height,h_iso21,rho_iso21)
		T = T_t2

	elif h_iso21 > height > h_iso12:
		T,p,rho = grad(p_iso12,rho_iso12,T_t1,a2,h_iso12,height)

	elif h_iso12 >= height >= h_iso11:
		p,rho = iso(p_iso11,T_t1,height,h_iso11,rho_iso11)
		T = T_t1

	elif h_iso11 > height > 0:
		T,p,rho = grad(p_sl,rho_sl,T_sl,a1,0,height)

	# NEGATIVE GEOMETRIC ALTITUDES ARE NOT SUPPORTED (Below sea level)

	# save each value of pressure, density and temperature on every iteration
	P.append(p)
	RHO.append(rho)
	t.append(T)	

	# Show the data for each altitude in the console	
	print("Temperature: " + str(round(T,4)) + " K")
	print("Pressure: " + str(round(p,4)) + " N/m2")
	print("Density: " + str(round(rho,6)) + " kg/m3\n\n--------------------------------")

# Plots for pressure, density and temperature vs. altitude
print(P,RHO,t)
plt.figure(1); plt.subplots_adjust(wspace=0, hspace=0.8)
plt.subplot(311); plt.plot(P,divide(h,1000),'r-', linewidth=0.8); plt.xlabel('Pressure (Pa)'); plt.ylabel('Altitude (km)'); plt.grid(True)
plt.subplot(312); plt.plot(RHO,divide(h,1000),'b-', linewidth=0.8); plt.xlabel('Density (kg/m3)'); plt.ylabel('Altitude (km)'); plt.grid(True)
plt.subplot(313); plt.plot(t,divide(h,1000),'g-', linewidth=0.8); plt.xlabel('Temperature (K)'); plt.ylabel('Altitude (km)'); plt.grid(True)
plt.show()