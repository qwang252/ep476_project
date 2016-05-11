""" NACA4_af class represents specific NACA airfoil shapes in 2D space"""
import numpy as np
import math
""" 4_symmetrical represents 4 digit unsymmetrical NACA airfoil""" 
class unsymmetrical(object):
  def __init__(self,name,is_FiniteTE,npanels=30):
    self.iaf = name
    self.is_FiniteTE = is_FiniteTE
    self.npanels = npanels    
  def indices_generate(self):
    # iaf is a string input
    m=int(self.iaf[0])/100
    p=int(self.iaf[1])/100
    t=int(self.iaf[2:4])/100
    
    a0= 0.2969
    a1=-0.1260
    a2=-0.3516
    a3= 0.2843
    
    if self.is_FiniteTE==1:
       a4=-0.1015
    else:
       a4=-0.1036
    
    beta=np.linspace(0,math.pi,self.npanels+1)
    x=0.5*(1-np.cos(beta))
    # x is a 2D column vector
    x=np.array([x]).T;
    # yt is column vector 
    yt=(t/0.2)*(a0*np.sqrt(x)+a1*x+a2*x**2+a3*x**3+a4*x**4)
    # xc1 is a 2D array  
    xc1=x[np.nonzero(x<=p)[0]]
    # xc2 is a 2D array 
    xc2=x[np.nonzero(x>p)[0]]
    # xc is 2D column vector
    xc=np.concatenate((xc1,xc2),axis=0)
    
    if p==0:
       xu=x
       yu=yt
       xl=x
       yl=-yt
       zc=np.zeros(np.shape(xc))
    else:
       yc1=(m/p**2)*(2*p*xc1-xc1**2) 
       yc2=(m/(1-p)**2)*((1-2*p)+2*p*xc2-xc2**2)
       zc=np.concatenate((yc1,yc2),axis=0)

       dyc1_dx =(m/p**2)*(2*p-2*xc1) 
       dyc2_dx =(m/(1-p)**2)*(2*p-2*xc2) 
       dyc_dx = np.concatenate((dyc1_dx,dyc2_dx),axis=0)       
       theta=np.arctan(dyc_dx)
       
       xu=x-yt*np.sin(theta) 
       yu=zc+yt*np.cos(theta) 
       
       xl=x+yt*np.sin(theta)
       yl=zc-yt*np.cos(theta)
       
       afx=np.concatenate((np.flipud(xu),xl[1::]),axis=0)
       afz=np.concatenate((np.flipud(yu),yl[1::]),axis=0)
       indices=np.array([[afx[i,0],afz[i,0]] for i in range(0,len(afx))])
       return indices
  def midpts(self,indices):
      n=len(indices)
      midpts=np.array([[(indices[0+i]+indices[-n+1+i])/2][0] for i in range(0,len(indices))])
      return midpts   
  def write_in_file(self,indices):
      # write data into file
       afxx=[]
       afzz=[]
       for i in range(0,len(indices)):
           afxx.append(indices[i,0])
           afzz.append(indices[i,1])
       
       Write_File_Name = 'NACA_'+self.iaf +'_airfoil.'+'dat'
       target = open(Write_File_Name, 'w')
       target.write('x indices'+'        '+'z indices'+'\n')
       for i in range(0,len(afxx)):
           target.write(str(afxx[i])+'    '+str(afzz[i]))
           target.write('\n')

      

