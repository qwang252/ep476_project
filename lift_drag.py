import numpy as np
import math
from flat_panel_velocity import flat_panel_velocity as fpv
import naca4_generator as naca
from pressure_calculation import pressure_cal as pc

def source_panel_LD(iaf,is_FiniteTE,npanels,rho_infinity,P_infinity,v_infinity,angle_attack):
    
    v_freestream=v_infinity*np.array([math.cos(angle_attack),math.sin(angle_attack),0])
    indices = normal_plotter(iaf,is_FiniteTE,npanels)[0]
    midpt = normal_plotter(iaf,is_FiniteTE,npanels)[1]
    n_hat = normal_plotter(iaf,is_FiniteTE,npanels)[2]
    t_hat = normal_plotter(iaf,is_FiniteTE,npanels)[3]    
    #normal velocity 
    N = normal_velocity(indices,midpt,n_hat)    
    print(n_hat.shape)
    n=len(indices)    
    #T_velo = np.zeros(shape=(npanels,npanels))
  
    B = np.dot(-v_freestream,n_hat.T)
    Ninv = np.linalg.inv(N-np.eye(n)/2)
    # lmda.shape = (1,61)
    lmda = np.dot(B,Ninv)
    #print(lmda.shape)
    #tangential velocity
    T = tangent_velocity(indices,midpt,t_hat,lmda)
   
    
    #tangential velocity for each panel

    # t_velo.shape = (1,61)
    t_velo = T.sum(axis = 0) + np.dot(v_freestream,t_hat.T)
    #pressure coefficient
    Cp = pressure_coeff(t_velo,v_infinity)
    length_panel = normal_plotter(iaf,is_FiniteTE,npanels)[4]
    
    F_net = net_force_cal(rho_infinity,v_infinity,P_infinity,Cp,length_panel,n_hat)
    
    LIFT = sum(np.dot(F_net,np.array([0,1,0]).T))
    DRAG = sum(np.dot(F_net,np.array([1,0,0]).T))
       
    return LIFT,DRAG
def normal_plotter(iaf,is_FiniteTE,npanels):
    
    indices = naca.naca4gen(iaf,is_FiniteTE,npanels)
    n=len(indices)
    indices = np.array([[indices[i,0],indices[i,1],0] for i in range(0,n)])
    
    #define normal and tangential directions for each panel
    indices_2 = np.concatenate((indices,np.array([indices[0]])),axis=0)
    
    length_panel = np.linalg.norm(np.diff(indices_2,axis=0),axis=1)
    t_hat = np.dot(np.diag(1/length_panel),np.diff(indices_2,axis=0))

    z_hat = np.array([0,0,1])
    z_hat = np.tile(z_hat,(n,1))
    print(z_hat)
    n_hat = np.cross(z_hat,t_hat)
    #define midpoints for each panel
    midpt = naca.midpts(indices)
    return indices,midpt,n_hat,t_hat,length_panel
    
def unit_tangent_dir(pos1,pos2):
    t_hat=(pos2-pos1)/np.linalg.norm(pos2-pos1)
    return t_hat

def normal_velocity(indices,midpt,n_hat):
    n=len(indices)
    L=[]
    for j in range(0,n):
        induced_vel_at_j =[np.dot(fpv(indices[j],indices[-n+1+j],midpt[i],1),n_hat[i]) for i  in range(0,n)]
        L.append(induced_vel_at_j)
    A = np.array(L)
    np.fill_diagonal(A,0) 
    return A
    
def tangent_velocity(indices,midpt,t_hat,lmda):
    n=len(indices)
    L=[]
    for j in range(0,n):
        induced_vel_at_j =[np.dot(fpv(indices[j],indices[-n+1+j],midpt[i],lmda[i]),t_hat[i]) for i  in range(0,n)]
        L.append(induced_vel_at_j)
    T = np.array(L)
    np.fill_diagonal(T,0)
    return T

def pressure_coeff(t_velo,v_infinity):
    Cp = 1-t_velo**2/v_infinity**2
    return Cp

def net_force_cal(rho_infinity,v_infinity,P_infinity,Cp,length_panel,n_hat):
    P = 0.5*Cp*rho_infinity*v_infinity**2+P_infinity*np.ones(shape=Cp.shape)
    P = np.diag(P)
    length = np.diag(length_panel)
    # n*n 
    Fp_mag = np.dot(length,P)
    F_net_panel = np.dot(Fp_mag,n_hat)
    return F_net_panel
    

    
