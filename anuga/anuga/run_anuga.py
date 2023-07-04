"""
Runs ANUGA through a BMI
"""

from __future__ import print_function

import sys
import numpy as np
from anuga_bmi import BmiAnuga


if __name__ == '__main__':

    sww = BmiAnuga()
    sww.initialize('anuga.yaml')
    

    for time in np.arange(0., 1000., 10):
    
        sww.update_until(time)
