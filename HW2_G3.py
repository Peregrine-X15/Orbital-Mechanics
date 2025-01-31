import numpy as np
import math
import cmath
import matplotlib.pyplot as plt

f_i=math.radians(-2.644) #rads
# omega=math.radians(329.0919)
# Omega=math.radians(351.2195)
ecc= 0.0253733
# i=math.radians(55.0856)
M_A_0=math.radians(209.8961)
n_revs=2.00565164156740
E_0=math.radians(3.651)

#n_rads is rads/day
#n_rads is rads/day
n_rads=n_revs*2*math.pi
n_radsPsec=n_rads/86400
print('n_radspSec=',n_radsPsec,'\n')
#assumed mass in kg
m_sat=1000
m_Earth=5.97219*(10**24)
G=6.67*(10**(-11))

mu=G*(m_sat+m_Earth) #mu in units of m^3/(kg*s^2)
# print('mu=',mu,'\n')
#from battin Lec 1
a=(mu/n_radsPsec**2)**(1/3)
param=a*(1-ecc**2)
r_0_I=[-26699.61973114e3,4792.5741857e3,946.61802299e3]# in m
v_0_I=[93096324.33000572*(10/36),19245110.87275932*(10/36),47609547.96201759*(10/36)]  # in m/s

del_t=12*3600
sigma_0=np.dot(r_0_I,v_0_I)/math.sqrt(mu)
print('sigma_0=',sigma_0,'\n')
M_A_t=M_A_0+n_radsPsec*del_t
print('M_A_t= ',M_A_t,'\n')

#newton's method
itr=51
E_t=np.zeros(itr)
E_t[0]=M_A_t
# Res_cond=0
Res_cond=1*10**-6
R2=1

for j in range(1, itr):
    E_t[j]=M_A_t+ecc*math.sin(E_t[j-1])
    R2=M_A_t*((ecc**j)/(1-ecc))-(E_t[j]-E_t[j-1])
    if(R2>=Res_cond):
        print('E_t=', E_t)
        print('R2=', R2)
        break
    else:
        continue
r_0_I_mag=math.sqrt(r_0_I[0]**2+r_0_I[1]**2+r_0_I[2]**2)
# print('r_0_I_mag=',r_0_I_mag)
print('\n')
E_t=E_t[1]
E_hat=E_t-E_0

r_t_mag=a+(r_0_I_mag-a)*math.cos(E_hat)+sigma_0*math.sqrt(a)*math.sin(E_0)
print('r_t_mag=', r_t_mag,'\n')
F=1-a/r_0_I_mag*(1-math.cos(E_hat))
G=del_t+math.sqrt(a**3/mu)*(math.sin(E_hat)-E_hat)
F_dot=-math.sqrt(mu*a)/(r_t_mag*r_0_I_mag)*math.sin(E_hat)
G_dot=1-a/r_t_mag*(1-math.cos(E_hat))


r_0_I=np.array([[-26699.61973114e3,4792.5741857e3,946.61802299e3]])# in m
v_0_I=np.array([[93096324.33000572*(10/36),19245110.87275932*(10/36),47609547.96201759*(10/36)]])  # in m/s
r_t=F*r_0_I+G*v_0_I
v_t=F_dot*r_0_I+G_dot*v_0_I
print('r_t=', r_t/1000,'km','\n')
print('v_t=', v_t*(3600/1000),'km/hr','\n')