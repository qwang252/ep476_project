# this funciton is to generate the NACA 4 digit airfoil coordinates with 
# specified number of panels (default number of panels is 30)
import numpy as np
import math
def naca4gen(iaf,is_FiniteTE,npanels=30):
    # iaf is a string input
    m=int(iaf[0])/100
    p=int(iaf[1])/100
    t=int(iaf[2:4])/100
    
    a0= 0.2969
    a1=-0.1260
    a2=-0.3516
    a3= 0.2843
    
    if is_FiniteTE==1:
       a4=-0.1015
    else:
       a4=-0.1036
    
    beta=np.linspace(0,math.pi,npanels+1)
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
       
       name=['NACA']
       afx=np.concatenate((np.flipud(xu),xl[1::]),axis=0)
       afz=np.concatenate((np.flipud(yu),yl[1::]),axis=0)

       # upper surface indices
       xU=afx[0:np.min(np.nonzero(afx==min(afx)))+1]
       #return np.min(np.nonzero(afx==min(afx)))
       zU=afz[0:np.min(np.nonzero(afx==np.min(afx)))+1]
       xL=afx[np.min(np.nonzero(afx==np.min(afx))):len(afx)+1]
       zl=afz[np.min(np.nonzero(afx==np.min(afx))):len(afx)+1]
       
       afxC=xc 
       afzC=zc
       
       lecirFactor=0.8
       afrLE=0.5*(a0*t/0.2)**2

       le_offs=0.5/100
       dyc_dx_le=(m/p**2)*(2*p-2*le_offs)
       #plt.plot(xU,zU)
       #xLEcenter=afrLE*math.cos(theta_le)
       #yLEcenter=afrLE*math.sin(theta_le)
       
       # write data into file
       FileHeader='NACA '+iaf+' airfoil'+'.dat'
       
       afxx=[]
       afzz=[]
       for i in range(0,len(afx)):
           afxx.append(afx[i,0])      
           afzz.append(afz[i,0])
       afxx=np.array(afxx)
       afzz=np.array(afzz)
       Write_File_Name = 'NACA_'+iaf +'_airfoil.'+'dat'
       target = open(Write_File_Name, 'w')
       target.write('x indices'+'        '+'z indices'+'\n')
       for i in range(0,len(afxx)):
           target.write(str(afxx[i])+'    '+str(afzz[i]))
           target.write('\n')
           
       return afxx, afzz     
       
            
       
