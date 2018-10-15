import numpy as np
import matplotlib.pyplot as plt

a=float(input('Enter a number: '))
b=float(input('Enter another number: '))
c=float(input('Enter another number: '))
d=float(input('Enter another number: '))
def func(x):
    return a*x**3+b*x**2+c*x+d

n=201
xvals=np.linspace(-10,10,n)
yvals=[]
for i in xvals:
    yvals.append(func(i))
    
plt.plot(xvals,yvals, 'r--')
plt.axhline(0,color='black')
plt.axvline(0,color='black')
plt.show
