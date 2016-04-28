import numpy as np
import math
from flat_panel_velocity import flat_panel_velocity as fpv
import naca4_generator as naca
#non-lifting flow around a cylinder of radius R with free stream
def source_panel(iaf,is_FiniteTE,npanels,v_infinity):
    
    v_freestream=v_infinity*np.array([[1,0,0]])

    #index of end poitns(xm,ym)
    indices = naca.naca4gen(iaf,is_FiniteTE,npanels)
    n=len(indices)
    indices = np.array([[indices[i,0],indices[i,1],0] for i in range(0,n)])
    
    #define normal and tangential directions for each panel
    t_hat=np.array([[unit_tangent_dir(indices[i],indices[-n+1+i])][0] for i in range(0,n)])    
    
    
    z_hat = np.array([0,0,1])
    z_hat = np.tile(z_hat,(n,1))
    
    n_hat = np.cross(z_hat,t_hat)
    
    
    #define midpoints for each panel
    midpt = naca.midpts(indices)
    
    #normal velocity 
    N = normal_velocity(indices,midpt,n_hat)    
    
    
    #T_velo = np.zeros(shape=(npanels,npanels))
  
    B = np.dot(-v_freestream,n_hat.T)
    Ninv = np.linalg.inv(N-np.eye(n)/2)
    # lmda.shape = (1,61)
    lmda = np.dot(B,Ninv)
    
    #tangential velocity
    T = tangent_velocity(indices,midpt,t_hat,lmda)
    
    
    #tangential velocity for each panel

    # t_velo.shape = (1,61)
    t_velo = T.sum(axis = 0) + np.dot(v_freestream,t_hat.T)
    
    #pressure coefficient
    Cp = pressure_coeff(t_velo,v_infinity)
    
 

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
        induced_vel_at_j =[np.dot(fpv(indices[j],indices[-n+1+j],midpt[i],lmda[0,i]),t_hat[i]) for i  in range(0,n)]
        L.append(induced_vel_at_j)
    T = np.array(L)
    np.fill_diagonal(T,0)
    return T

def pressure_coeff(t_velo,v_infinity):
    Cp = 1-t_velo**2/v_infinity**2
    return Cp
