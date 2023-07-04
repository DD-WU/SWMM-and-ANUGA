"""
Runs ANUGA through a BMI
"""

from __future__ import print_function

import math
import sys
import pandas as pd
from scipy.spatial import Delaunay
import anuga
import numpy as np
import matplotlib.pyplot as plt
from anugaBMI import BmiAnuga
from anuga.operators.set_stage_operator import Set_stage_operator
if __name__ == '__main__':
    import time as time1

  # 记录开始时间

    sww = BmiAnuga()
    sww.initialize_ANUGA('jiangbei.yaml')
    time_start = time1.time()
    print(time_start)
    for time in np.arange(0., 14400., 10):
        sww.update_until(time)
    # function()   执行的程序
    time_end = time1.time()
    print(time_end)
    time_sum = time_end - time_start
    print(time_sum)