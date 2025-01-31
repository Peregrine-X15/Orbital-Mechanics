import numpy as np
import math
import cmath
import matplotlib.pyplot as plt

#given parameters
ecc=0.25
a=12000E3
Omega=math.radians(260)
omega=math.radians(130)
f_o=math.radians(90)
i=math.radians(30)

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