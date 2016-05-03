from nose.tools import assert_equal
from naca4_generator as naca
def test_airfoil_thickness():
    airfoil_4415 = naca.naca4gen('4415',1,40)
    X = airfoil_4415[:,0]
    Y = airfoil_4415[:,1]
    obs = Y.max()-Y.min()
    exp = 0.15
    assert_equal(obs,exp)

def test_airfoil_camber()
    airfoil_4415 = naca.naca4gen('4415',1,40)
    XU = airfoil_4415[:,0][0:40]
    XL = airfoil_4415[:,0][41:]
    YU = airfoil_4415[:,1][0:40]
    YL = airfoil_4415[:,1][41:]
    

