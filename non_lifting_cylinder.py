import numpy as np
import math

#non-lifting flow around a cylinder of radius R with free stream
def non_lifting_cylinder(R,npanels,v_infinity):

    v_freestream=v_infinity*np.array([1,0,0])

    #define the panels as well as end points using complex
    #Npanels represents the number of panels, and therefore defines
    #discricized shape of cylinder

    dangle=2*math.pi/npanels
    
    #arrays/matrices for storing value of each loop
    theta = np.zeros(npanels)
    x = np.zeros(shape=(npanels,npanels))
    y = np.zeros(shape=(npanels,npanels))
    c_p = np.zeros(npanels)
    t = np.zeros(npanels + 1)
    t_hat = np.zeros(shape=(3,npanels))
    midpt = np.zeros(shape=(3,npanels))
    for m in range(1,npanels):
        theta[m-1] = dangle*(m - 1)-dangle/2

    #index of end poitns(xm,ym)
    for m in range(1,npanels):
        z = R*np.exp(1j*theta)
        x[:,m-1] = z.real
        y[:,m-1] = z.imag*1j

    #define normal and tangential directions for each panel
    for m in range(1,npanels):
        if m<=npanels - 1:
            t_hat[:,m - 1] = np.array([x[1,m] - x[1,m-1],y[1,m] - y[1,m - 1],0])/math.sqrt((x[1,m] - x[1,m-1])*(x[1,m] - x[1,m-1]) + (y[1,m] - y[1,m-1])*(y[1,m] - y[1,m-1]))
        else:
            t_hat[:,m - 1] = np.array([x[1,0] - x[1,m-1],y[1,0] - y[1,m - 1],0])/math.sqrt((x[1,0] - x[1,m-1])*(x[1,0] - x[1,m-1]) + (y[1,0] - y[1,m-1])*(y[1,0] - y[1,m-1]))


    #good
    z_hat = np.array([0,0,1])
    for m in range(1,npanels):
        n_hat = np.cross(z_hat,t_hat)


    #define midpoints for each panel
    for m in range(1,npanels):
        if m<= 1:
            midpt[:,m - 1] = np.array([x[1,m] +x[1,m - 1],y[1,m] + y[1,m - 1],0])/2
        else:
            midpt[:,m - 1] = np.array([x[1,0] + x[1,m - 1],y[1,0] + y[1,m - 1],0])/2



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
 
non_lifting_cylinder(1,2,3) 
