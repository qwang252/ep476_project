from nose.tools import assert_equal
import NACA_airfoil as naca
import numpy as np
import source_panel as sp
import vortex_panel as vp
airfoil_4415 = naca.unsymmetrical('4415',1,30)
def test_cl_source():
    obs = sp.source_panel(airfoil_4415,1,1,1,0)[0]
    exp = 0
    assert_equal(obs,exp)

def test_cd_source():
    obs = sp.source_panel(airfoil_4415,1,1,1,0)[1]
    exp = 0
    assert_equal(obs,exp)

def test_cl_vortex():
    obs = vp.vortex_panel(airfoil_4415,1,1,1,0)[0]
    exp = 0.5
    assert_equal(obs,exp)

def test_cd_vortex():
    obs = vp.vortex_panel(airfoil_4415,1,1,1,0)[1]
    exp = 0
    assert_equal(obs,exp)

