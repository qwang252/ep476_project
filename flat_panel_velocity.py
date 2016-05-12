
"""flat_panel_velocity is to calculate the velocity at desired position induced by a flat 
   source panel from pos1 to pos2 of strength lambda_1
   or to calculate the velocity at desired postion induced by a flat vortex panel 
   from pos1 to pos2 of strength gamma_1
...moduleauthor: Qingquan Wang<qwang252@wisc.edu> Yang Lou<lou9@wisc.edu>"""
import numpy as np
import math
def source_panel_velocity(pos1,pos2,desd_pos,lambda_1):
  # set position 1 in numpy array 
    pos1=np.array([pos1[0],pos1[1],0])
  # set position 2 in numpy array
    pos2=np.array([pos2[0],pos2[1],0])
  # set desired position in numpy array
    desd_pos=np.array([desd_pos[0],desd_pos[1],0])
    r1=desd_pos-pos1
    r2=desd_pos-pos2   
    r1_mag=np.linalg.norm(r1)
    r2_mag=np.linalg.norm(r2)
    p1p2=pos2-pos1
    p1p2_mag=np.linalg.norm(p1p2)
    t_hat=1/p1p2_mag*(p1p2)
    z_hat=np.array([0,0,1])
    r1r2_cross=np.cross(r1,r2)
    alpha=math.atan2(np.dot(r1r2_cross,z_hat),np.dot(r1,r2))
    n_hat=np.cross(z_hat,t_hat)
    
   # governing velocity equation for source panel method
    source_velocity=lambda_1/2/math.pi*(np.log(r1_mag/r2_mag)*t_hat+alpha*n_hat)
    return source_velocity 

def vortex_panel_velocity(pos1,pos2,desd_pos,gamma_1):
   # set position 1 in numpy array 
    pos1=np.array([pos1[0],pos1[1],0])
   # set position 2 in numpy array
    pos2=np.array([pos2[0],pos2[1],0])
   # set desired position in numpy array
    desd_pos=np.array([desd_pos[0],desd_pos[1],0])
    r1=desd_pos-pos1
    r2=desd_pos-pos2
    r1_mag=np.linalg.norm(r1)
    r2_mag=np.linalg.norm(r2)
    p1p2=pos2-pos1
    p1p2_mag=np.linalg.norm(p1p2)
    t_hat=1/p1p2_mag*(p1p2)
    z_hat=np.array([0,0,1])
    r1r2_cross=np.cross(r1,r2)
    alpha=math.atan2(np.dot(r1r2_cross,z_hat),np.dot(r1,r2))
    n_hat=np.cross(z_hat,t_hat)

   # governing velocity equation for vortex panel method
    vortex_velocity=gamma_1/2/math.pi*(-np.log(r1_mag/r2_mag)*n_hat+alpha*t_hat)
    return vortex_velocity
