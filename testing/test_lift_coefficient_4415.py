import lift_drag as ld
import numpy as np
def test_cl(ni = 5, nf = 40):
    cl_list = []
    for i in range(ni,nf):
        cl=abs(ld.lift_coefficient('4415',1,i,1,101325,1,0))
        cl_list.append(cl)
    cl_array=np.array(cl_list)
    return cl_array
    
