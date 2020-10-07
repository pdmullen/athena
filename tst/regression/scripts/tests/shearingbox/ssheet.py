# Regression test of shearing box and orbital advection with 2d hydro shwave.

# Modules
import logging
import math, cmath
import numpy as np
import scripts.utils.athena as athena
import sys
sys.path.insert(0, '../../vis/python')
import athena_read  # noqa
from mpmath import *
athena_read.check_nan_flag = True
logger = logging.getLogger('athena' + __name__[7:])  # set logger name based on module


# Prepare Athena++
def prepare(**kwargs):
    logger.debug('Running test ' + __name__)
    athena.configure(prob='ssheet',
                     eos='isothermal', **kwargs)
    athena.make()


# Run Athena++ w/wo Orbital Advection
def run(**kwargs):
    # w/o Orbital Advection
    arguments = [
        'job/problem_id=ssheet',
        'output1/file_type=hst', 'output1/dt=10.0',
        'output2/file_type=vtk', 'output2/variable=prim',
        'output2/dt=1000.0', 'time/cfl_number=0.4',
        'time/tlim=8000.0', 'time/nlim=2000',
        'time/xorder=2', 'time/integrator=vl2', 'time/ncycle_out=10',
        'mesh/nx1=64', 'mesh/x1min=-2.0', 'mesh/x1max=2.0',
        'mesh/ix1_bc=shear_periodic', 'mesh/ox1_bc=shear_periodic',
        'mesh/nx2=64', 'mesh/x2min=-2.0', 'mesh/x2max=2.0',
        'mesh/ix2_bc=periodic', 'mesh/ox2_bc=periodic',
        'mesh/nx3=1', 'mesh/x3min=-0.5', 'mesh/x3max=0.5',
        'hydro/iso_sound_speed=0.001', 'problem/ipert=1',
        'problem/amp=4.0e-4', 'problem/nwx=-4', 'problem/nwy=1', 'problem/nwz=0',
        'problem/Omega0=1.0e-3', 'problem/qshear=1.5', 'problem/shboxcoord=1',
        'problem/orbital_advection=false', 'time/ncycle_out=0']
    athena.run('hydro/athinput.ssheet', arguments)

    # w/  Orbital Advection
    arguments = [
        'job/problem_id=ssheet_oa',
        'output1/file_type=hst', 'output1/dt=10.0',
        'output2/file_type=vtk', 'output2/variable=prim',
        'output2/dt=1000.0', 'time/cfl_number=0.4',
        'time/tlim=8000.0', 'time/nlim=500',
        'time/xorder=2', 'time/integrator=vl2', 'time/ncycle_out=10',
        'mesh/nx1=64', 'mesh/x1min=-2.0', 'mesh/x1max=2.0',
        'mesh/ix1_bc=shear_periodic', 'mesh/ox1_bc=shear_periodic',
        'mesh/nx2=64', 'mesh/x2min=-2.0', 'mesh/x2max=2.0',
        'mesh/ix2_bc=periodic', 'mesh/ox2_bc=periodic',
        'mesh/nx3=1', 'mesh/x3min=-0.5', 'mesh/x3max=0.5',
        'hydro/iso_sound_speed=0.001', 'problem/ipert=1',
        'problem/amp=4.0e-4', 'problem/nwx=-4', 'problem/nwy=1', 'problem/nwz=0',
        'problem/Omega0=1.0e-3', 'problem/qshear=1.5', 'problem/shboxcoord=1',
        'problem/orbital_advection=true', 'time/ncycle_out=0']
    athena.run('hydro/athinput.ssheet', arguments)


# Analyze outputs
def analyze():
    # set parameters 
    ky     = 0.5*pi
    kx0    = -2.0*pi
    Omega0 = 0.001
    qshear = 1.5
    kappa2 = 2.0*(2.0-qshear)*Omega0**2.0
    cs     = 0.001
    constC = 0.5*(cs**2.0*ky**2.0+kappa2)/(qshear*Omega0*cs*ky)

    c1 = -1.82088e-07
    c2 = -8.20766e-08
    dvy0 = 1.0e-7

    # read results w/o Orbital Advection
    fname = 'bin/ssheet.hst'
    a    = athena_read.hst(fname)
    time1 = a['time']
    dvyc1 = a['dvyc']
    nf1   = len(time1)
    norm1 = 0.0
    for n in xrange(nf1):
        tau_    = qshear*Omega0*time1[n]+kx0/ky
        T_      = 1.0j * cmath.sqrt(2.0*cs*ky/(qshear*Omega0))*tau_
        exp_    = cmath.exp(-0.25j*T_*T_)
        fterm_  = exp_*hyp1f1(0.25-0.5j*constC, 0.5, 0.5j*T_*T_)
        first_  = fterm_.real
        exp_    = cmath.exp(0.25j*(pi-T_*T_))
        sterm_  = exp_*T_*hyp1f1(0.75-0.5j*constC, 1.5, 0.5j*T_*T_)
        second_ = sterm_.real
        advy = c1*first_+c2*second_
        norm1 += abs(dvyc1[n]-advy)/dvy0
    norm1 /= nf1

    # read results w/  Orbital Advection
    fname = 'bin/ssheet_oa.hst'
    b    = athena_read.hst(fname)
    time2 = b['time']
    dvyc2 = b['dvyc']
    nf2   = len(time2)
    norm2 = 0.0
    for n in xrange(nf2):
        tau_    = qshear*Omega0*time2[n]+kx0/ky
        T_      = 1.0j * cmath.sqrt(2.0*cs*ky/(qshear*Omega0))*tau_
        exp_    = cmath.exp(-0.25j*T_*T_)
        fterm_  = exp_*hyp1f1(0.25-0.5j*constC, 0.5, 0.5j*T_*T_)
        first_  = fterm_.real
        exp_    = cmath.exp(0.25j*(pi-T_*T_))
        sterm_  = exp_*T_*hyp1f1(0.75-0.5j*constC, 1.5, 0.5j*T_*T_)
        second_ = sterm_.real
        advy = c1*first_+c2*second_
        norm2 += abs(dvyc2[n]-advy)/dvy0
    norm2 /= nf2

    msg = '[ssheet]: L1 Norm of vy deviation {} = {}'
    logger.warning(msg.format('w/o Orbital Advection', norm1))
    logger.warning(msg.format('w/  Orbital Advection', norm2))
    flag = True
    if norm1 > 0.2:
        logger.warning('[SSHEET]: deviation is more than 20% w/o Orbital Advection')
        flag = False
    if norm2 > 0.2:
        logger.warning('[SSHEET]: deviation is more than 20% w/  Orbital Advection')
        flag = False

    return flag