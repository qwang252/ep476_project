from nose.tools import assert_equal
import NACA_airfoil as naca
import numpy as np
def test_airfoil_thickness():
    airfoil_4415 = naca.unsymmetrical('4415',1,30)
    X = airfoil_4415.indices_generate()[:,0]
    Y = airfoil_4415.indices_generate()[:,1]
    obs = Y.max()-Y.min()
    exp = 0.15
    assert_equal(obs,exp)

def test_airfoil_max_camber():
    airfoil_4415 = naca.unsymmetrical('4415',1,30)
    XU = airfoil_4415.indices_generate()[:,0][0:40]
    XL = airfoil_4415.indices_generate()[:,0][41:]
    YU = airfoil_4415.indices_generate()[:,1][0:40]
    YL = airfoil_4415.indices_generate()[:,1][41:]
    obs = max(YU - YL)
    exp = 0.04
    assert_equal(obs,exp)

def test_airfoil_location_maxcamber():
    airfoil_4415 = naca.unsymmetrical('4415',1,30)
    XU = airfoil_4415.indices_generate()[:,0][0:40]
    XL = airfoil_4415.indices_generate()[:,0][41:]
    YU = airfoil_4415.indices_generate()[:,1][0:40]
    YL = airfoil_4415.indices_generate()[:,1][41:]
    diff = YU - YL
    obs = np.argmax(diff)
    exp = 30*0.4
    assert_equal(obs,exp)


