import numpy as np
import math
import cmath
import matplotlib.pyplot as plt

omega=math.radians(329.0919)
Omega=0
i=0

ROI=np.array([[math.cos(omega)*math.cos(Omega)-math.sin(omega)*math.cos(i)*math.sin(Omega), math.cos(omega)*math.sin(Omega)+math.sin(omega)*math.cos(i)*math.cos(Omega), math.sin(omega)*math.sin(i)],
    [-(math.sin(omega)*math.cos(Omega)+math.cos(omega)*math.cos(i)*math.sin(Omega)), math.cos(omega)*math.cos(i)*math.cos(Omega)-math.sin(omega)*math.sin(Omega), math.cos(omega)*math.sin(i)],
    [math.sin(i)*math.sin(Omega), -math.sin(i)*math.cos(Omega), math.cos(i)]])

print('ROI',ROI)
ROI_T=ROI.transpose()