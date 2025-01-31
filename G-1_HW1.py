import numpy as np
import math
import cmath
import matplotlib.pyplot as plt
# from Planet3D import Planet, generate_solar_system
#from orbital import earth, KeplerianElements, Maneuver


#pi=3.141592653589793
#a=np.linspace(0,50,70)
#print('a=' , a)

#Given values from TLE
ecc= 0.0253733
inclin=math.radians(55.0856)
omega=math.radians(329.0919)
Omega=math.radians(351.2195)
i=math.radians(55.0856)
M_A=math.radians(209.8961)
print('M_A=',M_A,'\n')
n_revs=2.00565164156740
#n_rads is rads/day
n_rads=n_revs*2*math.pi
print('n_rads=',n_rads,'rads/day','\n')
print('n_rads=',n_rads/86400,'rads/sec','\n')
n_radsPsec=n_rads/86400
#assumed mass in kg
m_sat=1000
m_Earth=5.97219*(10**24)
G=6.67*(10**(-11))

mu=G*(m_sat+m_Earth) #mu in units of m^3/(kg*s^2)
print('mu=',mu,'\n')
#from battin Lec 1
a=(mu/n_radsPsec**2)**(1/3)
print('a= ', a/1000, 'Km','\n')

# Period_rads=2*math.pi/n_rads
# Period_revs=2*math.pi/n_revs
# print('Period =',Period_rads, 'days per 2pi rads')
# print('Period =',Period_revs, 'days per revolution')

P=2*math.pi/math.sqrt((mu/a**3))
print('P= ', P, 'seconds')
print('P= ', P/3600, 'hours','\n')
#print('P= ', 2*math.pi/n_radsPsec, 'hours')
c=mu*ecc
print('c=',c)

# ROI=[[math.cos(omega)*math.cos(Omega)-math.sin(omega)*math.cos(i)*math.sin(Omega), math.cos(omega)*math.sin(Omega)+math.sin(omega)*math.cos(i)*math.cos(Omega), math.sin(omega)*math.sin(i)],
#     [-(math.sin(omega)*math.cos(Omega)+math.cos(omega)*math.cos(i)*math.sin(Omega)), math.cos(omega)*math.cos(i)*math.cos(Omega)-math.sin(omega)*math.sin(Omega), math.cos(omega)*math.sin(i)],
#     [math.sin(i)*math.sin(Omega), -math.sin(i)*math.cos(Omega), math.cos(i)]]
ROI=np.array([[math.cos(omega)*math.cos(Omega)-math.sin(omega)*math.cos(i)*math.sin(Omega), math.cos(omega)*math.sin(Omega)+math.sin(omega)*math.cos(i)*math.cos(Omega), math.sin(omega)*math.sin(i)],
    [-(math.sin(omega)*math.cos(Omega)+math.cos(omega)*math.cos(i)*math.sin(Omega)), math.cos(omega)*math.cos(i)*math.cos(Omega)-math.sin(omega)*math.sin(Omega), math.cos(omega)*math.sin(i)],
    [math.sin(i)*math.sin(Omega), -math.sin(i)*math.cos(Omega), math.cos(i)]])

ROI_T=ROI.transpose()
print('ROI= ', ROI,'\n')
itr=51
E=np.zeros(itr)
Res_cond=0
R2=1
#while R2>=Res_cond:
for j in range(1, itr):
    E[j]=M_A+ecc*math.sin(E[j-1])
    R2=M_A*((ecc**j)/(1-ecc))-(E[j]-E[j-1])
    if(R2>=Res_cond):
        print('E=', E)
        print('R2=', R2)
        break
    else:
        continue
E=E[2]
#print(E)
f=2*math.atan((math.sqrt(1+ecc)/math.sqrt(1-ecc))*math.tan(0.5*E))
print('f=', f,'\n')

param=a*(1-ecc**2)
print('p(paramter)=',param/1000,'Km','\n')

#part d
r_mag=param/(1+ecc*math.cos(f))
rox=r_mag*math.cos(f)
roy=r_mag*math.sin(f)
print('rox=',rox,'\n')
print('roy=',roy,'\n')
r_Omag_matrix=np.array([[rox],[roy],[0]])

#r_Imag_matrix=ROI_T*r_Omag_matrix
r_Imag_matrix=np.matmul(ROI_T,r_Omag_matrix)/1000
print('ROI=',ROI,'\n')
print('ROI_T', ROI_T,'\n')
print('r_Imag_matrix=',r_Imag_matrix,'\n','in Kms','\n')

E_dot=n_radsPsec/(1-ecc*math.cos(E))

v_O=np.array([[-a*math.sin(E)],[a*math.sqrt(1-ecc**2)],[0]])
print('v_0=',v_O)
v_I=np.matmul(ROI_T,v_O)
#print('v_I=',v_I,'m/s')
print('v_I=',v_I*(3600/1000),'km/hr','\n')

#part e
r_a=a*(1-ecc)
r_p=a*(1+ecc)

print('r_a=',r_a/1000,'Km','\n')
print('r_p=',r_p/1000,'Km','\n')

v_a=math.sqrt(2*(-mu/(2*a)+mu/r_a))
v_p=math.sqrt(2*(-mu/(2*a)+mu/r_p))

print('v_a=',v_a,'m/s','\n')
print('v_p=',v_p,'m/s','\n')

#part f
r_pf=20E3 #m
v_pf=math.sqrt(2*(-mu/(2*a)+mu/r_pf))
print('v_pf=',v_pf,'m/s','\n')

v_f_esc=math.sqrt(2*mu/r_pf)
print('v_f_esc=',v_f_esc,'m/s','\n')

del_v_req=v_f_esc-v_pf
print('Δv_req=',del_v_req,'m/s','\n')

#part g
#f= -2.6446879385485302  rads

#part h
f_p=0
r_op=np.array([[r_p*math.cos(f_p)],[r_p*math.sin(f_p)],[0]])
# print('r_op=',r_op,'m','\n')
r_Ip=np.matmul(ROI_T,r_op)
print('r_Ip=',r_Ip,'m','\n')
#part i
f_i=math.pi/2
r_pi=param/(1+ecc*math.cos(f_i))
r_o_pi_vec=np.array([[r_pi*math.cos(f_i)],[r_pi*math.sin(f_i)],[0]])
r_i_pi=np.matmul(ROI_T,r_o_pi_vec)
print('r_i_pi=',r_i_pi,'m','\n')
#part j

a_g_p=mu/r_p**2
a_g_a=mu/r_a**2

print('acceleration due to gravity at periapsis',a_g_p,'m/s^2','\n')
print('acceleration due to gravity at apoapsis',a_g_a,'m/s^2','\n')

######################################################################################
#G-2
#part c
omega_bar=Omega+omega
print('omega_bar=',omega_bar,'rads')
nu=omega_bar+M_A
print('nu',nu,'rads')
print('a=',a/1000,'km','\n')
h_ii=ecc*math.sin(nu)
print('h_ii=',h_ii,)
k_ii=ecc*math.cos(nu)
print('k_ii=',k_ii)
p_ii=math.tan(i/2)*math.sin(Omega)
print('p_ii=',p_ii)
q_ii=math.tan(i/2)*math.cos(Omega)
print('q_ii',q_ii)

