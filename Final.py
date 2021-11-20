# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 16:29:42 2017

@author: Tom
"""

from physics import *
import pylab

#   number of electrons per nucleon
Y_e = 0.5

#   mass of a proton, in kg
m_p = 1.6726231*10**(-27.0)

#   mass of an electron, in kg
m_e = 9.10938291*10**(-31.0)

#   Planck constant, in m^2*kg/s
h = 6.62606957*10**(-34.0)

#   gravitational constant, in m^3 kg^-1 s^-2
G = 6.67384*10**(-11.0)

#   speed of light, in m/s
c = 2.99792458*10**(8.0)

#   core density of the sun, in kg/m^3
rho_sun = 1.622*10**(5.0)

#   mass of the sun, in kg
m_sun = 1.98855*10**(30.0)

#   radius of the sun, in m
r_sun = 6.955*10**(8.0)

#   radius of the earth, in m
r_earth = 6.3781*10**(6.0)

#   define rho_0^-1/3
H = ((3.0*Y_e)/(8.0*pi*m_p))**(1.0/3.0)*(h/(m_e*c))

#   receive rho (core density) input
rho_input = input("Enter core density value: ")

#   define gamma factor
def Y(rho_input, y, H):
    result = (H**(2.0)*rho_input**(2.0/3.0)*abs(y[1])**(2.0/3.0))/(3*(1+H**(2.0)*rho_input**(2.0/3.0)*abs(y[1])**(2.0/3.0))**(1.0/2.0))
    return result

#   define scale factors
c_m = (4*pi*r_sun**(3.0)*rho_input)/m_sun
c_rho = (-1*G*m_p*m_sun)/(Y_e*m_e*c**(2.0)*r_sun)
    
#----------------------------------------------------------------------------

def derivs(n, r, y):
	"The function DERIVS calculates y' from x and y"
	dy=[0 for i in range(0,n+1)]
	if r==0.0: dy[1]=0.0
	else:      dy[1]=(c_rho*y[2]*abs(y[1]))/(r**(2.0)*Y(rho_input, y, H))

	dy[2]=c_m*r**(2.0)*y[1]
	return dy
#----------------------------------------------------------------------------

def runkut(n, r, y, h):
	"Advances the solution of diff eqn defined by derivs from x to x+h"
	y0=y[:]
	k1=derivs(n, r, y)
	for i in range(1,n+1): y[i]=y0[i]+0.5*h*k1[i]
	k2=derivs(n, r+0.5*h, y)
	for i in range(1,n+1): y[i]=y0[i]+h*(0.2071067811*k1[i]+0.2928932188*k2[i])
	k3=derivs(n, r+0.5*h, y)
	for i in range(1,n+1): y[i]=y0[i]-h*(0.7071067811*k2[i]-1.7071067811*k3[i])
	k4=derivs(n, r+h, y)
	for i in range(1,n+1):
		a=k1[i]+0.5857864376*k2[i]+3.4142135623*k3[i]+k4[i]
		y[i]=y0[i]+0.16666666667*h*a
	
	r+=h
	return (r,y)
#---------------------------------------------------------------------------- 

r=0.0; y=[0, rho_input/rho_sun, 0]

N=10000000

w=[0 for j in range(0,N)]
z=[0 for j in range(0,N)]
q=[0 for j in range(0,N)]

j=0
while y[1]>0:
      (r,y) = runkut(2, r, y, 1.0/N)
      print (r, y[1], y[2])
      w[j]=y[1]
      z[j]=y[2]
      q[j]=r
      j=j+1
      
pylab.plot(q[:j],z[:j])                     #doesn't plot zero values
pylab.xlabel('Radius (Solar Units)')
pylab.ylabel('Mass (Solar Units)')
pylab.show()

pylab.plot(q[:j],w[:j])                     #doesn't plot zero values
pylab.xlabel('Radius (Solar Units)')
pylab.ylabel('Density (Solar Units)')
pylab.ylim(0)
pylab.show()