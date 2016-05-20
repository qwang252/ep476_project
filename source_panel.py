"""source_panel represents source panel method in aerodynamics
...moduleauthor::Qingquan Wang<qwang252@wisc.edu>, Yang Lou<lou9@wisc.edu>"""
from flat_panel_velocity import source_panel_velocity as fspv
import truncation_plotter as tp
import numpy as np
import math
# source_panel function returns cl and cd for specific airfoil
def source_panel(af,rho_infinity,P_infinity,v_infinity,angle_attack):
   # define free stream velocity with v_infinity and angle_attack
    v_freestream=v_infinity*np.array([math.cos(angle_attack),math.sin(angle_attack),0])
   # set indices,mid points,unit normal vector and unit tangential vector  
   # of each panel for specific airfoil
    indices = tp.truncation_plotter(af)[0]
    midpt = tp.truncation_plotter(af)[1]
    n_hat = tp.truncation_plotter(af)[2]
    t_hat = tp.truncation_plotter(af)[3]
   # normal velocity at the midpt for each panel 
    N = normal_velocity(indices,midpt,n_hat)

    n=len(indices)

    B = np.dot(-v_freestream,n_hat.T)
    Ninv = np.linalg.inv(N-np.eye(n)/2)
    lmda = np.dot(B,Ninv)
   # tangential velocity at the midpt for each panel
    T = tangent_velocity(indices,midpt,t_hat,lmda)


    t_velo = T.sum(axis = 0) + np.dot(v_freestream,t_hat.T)
   # pressure coefficient at the midpt for each panel
    Cp = pressure_coeff(t_velo,v_infinity)
    length_panel = truncation_plotter(af)[4]

    F_net = net_force_cal(rho_infinity,v_infinity,P_infinity,Cp,length_panel,n_hat)
   # net lift and drag force
    LIFT = sum(np.dot(F_net,np.array([0,1,0]).T))
    DRAG = sum(np.dot(F_net,np.array([1,0,0]).T))
    mean_y = np.mean(indices[:,1])
   # area of specific airfoil
    A_af = mean_y
    Cl = LIFT/rho_infinity/v_infinity**2*2
    Cd = DRAG/rho_infinity/v_infinity**2/A_af*2/mean_y
    return Cl,Cd

# calculate normal velocity at midpt of each panel
def normal_velocity(indices,midpt,n_hat):
    n=len(indices)
    L=[]
    for j in range(0,n):
         induced_vel_at_j =[np.dot(fspv(indices[j],indices[-n+1+j],midpt[i],1),n_hat[i]) for i in range(0,n)] 
         L.append(induced_vel_at_j)
    A=np.array(L)
    np.fill_diagonal(A,0)
    return A

# calculate tangential velocity at midpt of each panel
def tangent_velocity(indices,midpt,t_hat,lmda):
    n=len(indices)
    L=[]
    for j in range(0,n):
        induced_vel_at_j =[np.dot(fspv(indices[j],indices[-n+1+j],midpt[i],lmda[i]),t_hat[i]) for i in range(0,n)]
        L.append(induced_vel_at_j)
    T = np.array(L)
    np.fill_diagonal(T,0)
    return T

# calculate pressure coefficient at midpt of each panel
def pressure_coeff(t_velo,v_infinity):
    Cp = 1-t_velo**2/v_infinity**2
    return Cp

# integrate pressure force on each panel to calculate net force on the whole airfoil
def net_force_cal(rho_infinity,v_infinity,P_infinity,Cp,length_panel,n_hat):
    P = 0.5*Cp*rho_infinity*v_infinity**2+P_infinity*np.ones(shape=Cp.shape)
    P = np.diag(P)
    length = np.diag(length_panel)
    # n*n 
    Fp_mag = np.dot(length,P)
    F_net_panel = np.dot(Fp_mag,n_hat)
    return F_net_panel
