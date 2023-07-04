
from __future__ import print_function

import sys
import numpy as np
from anuga_bmi import BmiAnuga
import anuga

if __name__ == '__main__':

    sww = BmiAnuga()
    sww.initialize('jiangbei.yaml')
    center = (658538.453, 3550528.655)
    radius = 10.0
    region0 = anuga.Region(sww._anuga.domain, center=center, radius=radius)
    fixed_inflow = anuga.Inlet_operator(sww._anuga.domain, region0, Q=19.7)
    for time in np.arange(0., 1000., 10):
        fixed_inflow.set_Q(1000)
        sww.update_until(time)