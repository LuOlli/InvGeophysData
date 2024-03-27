import numpy as np
import matplotlib.pyplot as plt
x = np.array(
    [2.5,2,1,0]#x-coordinate data
)
y = np.array(
    [3.76,6.01,8.98,10.02]#y-coordinate data
)
one = np.ones((1,x.size)) #a vector of ones to represent the constant term
G = np.transpose(np.vstack((one,x,x**2))) #the observation matrix for fitting a parabola
GTG = np.matmul(np.transpose(G),G) #matrix G^TG
GTz = np.matmul(np.transpose(G),np.transpose(y)) #vector G^Tz
LS = np.matmul(np.linalg.inv(GTG),GTz) # the least-squares solution parameters
x0=np.roots([LS[2],LS[1],LS[0]])#the location of the thrower
print('Location of the thrower:',x0[1])
lp = np.linspace(x.min(),3.5,100) # a linspace vector for plotting the LS solution
m = LS[0]+LS[1]*lp+LS[2]*lp**2 # the least-squares solution
print('LS solution:',np.round(LS[2], 3),'*x^2+',np.round(LS[1], 3),'*x+',np.round(LS[0], 3))
plt.plot(x,y,'ks')
plt.plot(lp,m,'r')
plt.xlabel('x')
plt.ylabel('y') 
plt.ylim(0,10.5)
plt.savefig('BallThrow.png')
