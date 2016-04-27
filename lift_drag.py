import numpy as np
import math
import flat_panel_velocity
import naca4_generator as naca
#non-lifting flow around a cylinder of radius R with free stream
def non_lifting_cylinder(iaf,is_FiniteTE,npanels,v_infinity):
    
    v_freestream=v_infinity*np.array([1,0,0])

    #index of end poitns(xm,ym)
    indices = naca.naca4gen(iaf,is_FiniteTE,npanels)

    #define normal and tangential directions for each panel
    n=len(indices)
    t_hat=np.array([[unit_tangent_dir(indices[0+i],indices[-n+1+i])][0] for i in range(0,n)])    

    
    z_hat = np.array([0,0,1])
    z_hat = np.tile(z_hat,(n,1))
    
    n_hat = np.cross(z_hat,t_hat)

    print(n_hat.shape)
    #define midpoints for each panel
    indices_midpt = naca.midpts(indices)
    # everything until here is good
    #governing velocity equation
    A = np.zeros(shape=(npanels,npanels))
    T_velo = np.zeros(shape=(npanels,npanels))
    for i in range(1,npanels):
        for j in range(1,npanels):
            if i == j:
                if j < npanels:
                    A[j-1,i-1]=flat_panel_velocity(x[1,j-1],y[1,j-1],x[1,j],y[1,j],midpt[0,i-1],midpt[1,i-1],1)*n_hat[:,i-1]
                    
                else:
                    A[j-1,i-1]=flat_panel_velocity(x[1,j-1],y[1,j-1],x[1,0],y[1,0],midpt[0,i-1],midpt[1,i-1],1)*n_hat[:,i-1]



    lmda = np.zeros(npanels)
    B = -v_freestream*n_hat
    lmda = B*np.linalg.inv(A - np.eye(npanels)/2)

    #tangential velocity
    for i in range(1,npanels):
        for j in range(1,npanels):
            if i == j:
                if j < npanels:
                    T_velo[j-1,i-1]=flat_panel_velocity(x[1,j-1],y[1,j-1],x[1,j],y[1,j],midpt[0,i-1],midpt[1,i-1],lmda[0,j-1])*n_hat[:,i-1]
                    
                else:
                    T_velo[j-1,i-1]=flat_panel_velocity(x[1,j-1],y[1,j-1],x[1,0],y[1,0],midpt[0,i-1],midpt[1,i-1],lmda[0,j-1])*n_hat[:,i-1]



    #tangential velocity for each panel
    t_velo = T.sum(axis = 0) + v_freestream*t_hat

    #pressure coefficient
    for j in range(1,npanels + 1):
        if j <= npanels:
            c_p[j-1] = 1 - (t_velo[j-1]*t_velo[j-1])/(v_infinity*v_infinity)
        else:
            c_p[j-1] = 1 - (t_velo[j-1]*t_velo[0])/(v_infinity*v_infinity)
    
	#plot cp
    for i in range(1,npanels + 1):
        t[i-1] = dangle*(i - 2)
    
    tt = np.array([np.arange(0, 1, 1/math.pow(2, 6))])*2*math.pi
    c_p_an = 1 - 4*(np.sin(tt)*np.sin(tt))
 
#non_lifting_cylinder(1,2,3) 

def unit_tangent_dir(pos1,pos2):
    t_hat=(pos2-pos1)/np.linalg.norm(pos2-pos1)
    return t_hat

def gov_vel_eqn():
    pass 
