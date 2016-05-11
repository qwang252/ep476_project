
# this function is to calculate the velocity at desired position induced by a flat 
# source panel from pos1 to pos2 of strength lambda_1
# or to calculate the velocity at desired postion induced by a flat vortex panel 
# from pos1 to pos2 of strength gamma_1
import numpy as np
import math
def flat_sourece_panel_velocity(pos1,pos2,dest_pos,lambda_1): 
    pos1=np.array([pos1[0],pos1[1],0])
    pos2=np.array([pos2[0],pos2[1],0])
    dest_pos=np.array([dest_pos[0],dest_pos[1],0])
    r1=dest_pos-pos1
    r2=dest_pos-pos2   
    r1_mag=np.linalg.norm(r1)
    r2_mag=np.linalg.norm(r2)
    p1p2=pos2-pos1
    p1p2_mag=np.linalg.norm(p1p2)
    t_hat=1/p1p2_mag*(p1p2)
    z_hat=np.array([0,0,1])
    r1r2_cross=np.cross(r1,r2)
    alpha=math.atan2(np.dot(r1r2_cross,z_hat),np.dot(r1,r2))
    n_hat=np.cross(z_hat,t_hat)
    
    
    source_velocity=lambda_1/2/math.pi*(np.log(r1_mag/r2_mag)*t_hat+alpha*n_hat)
    return source_velocity 
def flat_vortex_panel_velocity(pos1,pos2,dest_pos,gamma_1):
    pos1=np.array([pos1[0],pos1[1],0])
    pos2=np.array([pos2[0],pos2[1],0])
    dest_pos=np.array([dest_pos[0],dest_pos[1],0])
    r1=dest_pos-pos1
    r2=dest_pos-pos2
    r1_mag=np.linalg.norm(r1)
    r2_mag=np.linalg.norm(r2)
    p1p2=pos2-pos1
    p1p2_mag=np.linalg.norm(p1p2)
    t_hat=1/p1p2_mag*(p1p2)
    z_hat=np.array([0,0,1])
    r1r2_cross=np.cross(r1,r2)
    alpha=math.atan2(np.dot(r1r2_cross,z_hat),np.dot(r1,r2))
    n_hat=np.cross(z_hat,t_hat)
    
    vortex_velocity = gamma_1/2/math.pi*(-n_hat*np.log(r1_mag/r2_mag)+t_hat*alpha)
    return vortex_velocity


