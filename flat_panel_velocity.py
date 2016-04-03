# this function is to calculate the velocity at (x,y) induced by a flat 
# source panel from (x1,y1) to (x2,y2) of strength lambda_1
import numpy as np
import math
def flat_panel_velocity(x1,y1,x2,y2,x,y,lambda_1):
    p1=np.array([x1,y1,0])
    p2=np.array([x2,y2,0])
    p=np.array([x,y,0])
    r1=p-p1
    r2=p-p2
    r1_mag=np.sqrt((r1*r1).sum())
    r2_mag=np.sqrt((r2*r2).sum())
    p1p2=p1-p2
    p1p2_mag=np.sqrt((p1p2*p1p2).sum())
    t_hat=1/p1p2_mag*(p2-p1)
    z_hat=np.array([0,0,1])
    alpha=math.atan2(np.dot(np.cross(r1,r2),z_hat),np.dot(r1,r2))
    n_hat=np.cross(z_hat,t_hat)
    velocity=lambda_1/2/math.pi*(np.log(r1_mag/r2_mag)*t_hat+alpha*n_hat)
    return velocity 

