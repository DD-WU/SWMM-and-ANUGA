#! /usr/bin/env python
"""Basic Model Interface implementation for anuga."""
import anuga
import numpy as np
import scipy
import yaml
from basic_modeling_interface import Bmi
from anuga.operators.set_stage_operator import Set_stage_operator
from anuga_solver import AnugaSolver
from anuga.operators.rate_operators import Rate_operator

class BmiAnuga(Bmi):

    def __init__(self):
        """Create a BmiAnuga model that is ready for initialization."""
        self._anuga = None
        self._time = 0.
        self._values = {}
        self._var_units = {}
        self._grids = {}
        self._grid_type = {}
        self.regions={}
        self.elevation={}
        self.rate_operator={}
        self.area_map={}
        self.rain={}
    """
    process interface
    """
    def initialize_ANUGA(self, filename='anuga.yaml'):
        """Initialize the ANUGA model.

        Parameters
        ----------
        filename : str, optional
            Path to name of input file.
        """

        with open(filename, 'r') as file_obj:
            params = yaml.load(file_obj,Loader=yaml.FullLoader)

        default_params = {'domain_type':'square',
                          'shape':(10.,5.),
                          'size':(10.,5.),
                          'friction':0.,
                          'boundary_filename':'',
                          'point_and_vertices':[],
                          'pipe_Index_Cord':'',
                          'rainfall': '',
                          'elevation_filename':'',
                          'elevation_profile':'shallow linear ramp',
                          'output_filename':'anuga_output',
                          'output_timestep':10,
                          "save_timestep":0,
                          'boundary_tags':{'left':[],
                                           'right':[],
                                           'top':[],
                                           'bottom':[]},
                          'boundary_conditions':{'left': 'Reflective',
                                       'right': ['Dirichlet', 5, 0, 0],
                                       'top': 'Reflective',
                                       'bottom': 'Reflective'},
                          'stored_quantities':{'stage':2,
                                               'xmomentum':2,
                                               'ymomentum':2,
                                               'elevation':1},
                          'maximum_triangle_area':10,
                          'regionPtArea_filename':'',
                          'initial_flow_depth': 0,
                          'interior_polygon_filename':'',
                          'interior_polygon_triangle_area': 0.0,
                          'toggle_sediment_transport':False,
                          'inflow_sediment_concentration': 0.0,
                          'initial_sediment_concentration': 0.0,
                          'toggle_vegetation_drag': False,
                          'vegetation_stem_diameter': 0.0,
                          'vegetation_stem_spacing': 0.0,
                          'Mannings_n_parameter': 0.0,
                          }

        for key,value in default_params.items():
            params[key] = params.get(key, default_params[key])

        for key,value in params.items():
            if (value is None) or (value == 'None'):
                params[key] = default_params[key]

        assert (set(params['boundary_conditions'].keys()) ==
                set(params['boundary_tags'].keys())), (
                "The boundary tag names don't match the boundary "
                "condition names. Check that the two dictionaries use "
                "the same boundary names.")

        self._anuga = AnugaSolver(params)
        self.create_rainfall()
        self.create_pipe_dict()
        self.create_elevation_dict()
        self.create_rate_dict()
        self.create_area_dict()
    def update_test(self):
        for t in self.domain.evolve(yieldstep = self._time_step, finaltime = self._time):
            print(self.domain.timestepping_statistics())
    def update_ANUGA(self):
        """Advance model by one time step."""
        self._time += self.get_time_step_ANUGA()
        self._anuga._time = self._time
        self._anuga.update()
    def update_frac_ANUGA(self, time_frac):
        """Update model by a fraction of a time step.

        Parameters
        ----------
        time_frac : float
            Fraction fo a time step.
        """
        time_step = self.get_time_step_ANUGA()
        self._anuga.time_step = time_frac * time_step
        self.update_ANUGA()
        self._anuga.time_step = time_step
    def update_until(self, then):
        """Update model until a particular time.
        Parameters----------then : float Time to run model until."""
        n_steps = (then - self.get_current_time_ANUGA()) / self.get_time_step_ANUGA()
        for _ in range(int(n_steps)):
            self.update_ANUGA()
        if (n_steps - int(n_steps)) > 0.0:
            self.update_frac_ANUGA(n_steps - int(n_steps))
    def finalize_ANUGA(self):
        """Finalize model."""
        self._anuga = None
    """
    data interface
    """
    def create_rainfall(self,delimiter=",",polygon=None):
        rain_timeseries = scipy.genfromtxt(
            self._anuga._rainfall, delimiter=',', skip_header=1)
        # Adjust starttime
        rain_timeseries[:, 0] = rain_timeseries[:, 0]

        # Convert units to m/s (from mm/hr)
        rain_timeseries[:, 1] = rain_timeseries[:, 1] / (3600. * 1000.)

        # Sanity check
        assert rain_timeseries[:, 1].min() >= 0., 'Negative rainfall input'

        # Make interpolation function and add to ANUGA as operator
        if rain_timeseries[:, 1].max() >= 0.:
            myrain = scipy.interpolate.interp1d(
                rain_timeseries[:, 0], rain_timeseries[:, 1],
                kind='cubic')
            anuga.operators.rate_operators.Rate_operator(
                self._anuga.domain, rate=myrain, polygon=polygon, label= self._anuga._rainfall)
        return

    def create_elevation_dict(self):
        for key,value in self.regions.items():
            self.elevation[key] = float(value.domain.quantities['elevation'].centroid_values[value.indices])
    def create_pipe_dict(self,radius=100):
        region_temp = open(self._anuga._pipe_Index_Cord)
        lines = region_temp.readlines()
        region_temp.close()
        t={}
        i=0
        for line in lines:
            fields = line.split(',')
            r=anuga.Region(self._anuga.domain,center=(float(fields[0]),float(fields[1])),radius=radius,capture=True)
            t[fields[2].strip()]= r
            i+=1
            if i==275:
                print()
            print(fields[2].strip())
        self.regions=t


    def create_rate_dict(self):
        # 为了方便监测rate_operator,因此这边注释掉，放到region中
        for key ,value in self.regions.items():
            rate = anuga.Rate_operator(self._anuga.domain, rate=0, indices=value.indices)
            self.rate_operator[key] = rate

    def create_area_dict(self):
        for key, value in self.regions.items():
            self.area_map[key]=float(value.areas[value.indices])
    def get_start_time_ANUGA(self):
        """Start time of model."""
        return 0.
    def get_end_time_ANUGA(self):
        """End time of model."""
        return np.finfo('d').max
    def get_current_time_ANUGA(self):
        """Current time of model."""
        return self._time
    def get_time_step_ANUGA(self):
        """Time step of model."""
        return self._anuga.time_step
    def set_stage_ANUGA(self,value=0,name=None,location="centroids",indices=None):
        if indices is None:
            indices=self.regions[name].indices
        self._anuga.domain.set_quantity("stage",location=location,indices=indices,numeric=self.elevation[name]+value)
    def set_rate_ANUGA(self, value=0, name=None, location="centroids", indices=None):
        self.rate_operator[name].set_rate(value)
    def get_stage_ANUGA(self,name=None,location="centroids",indices=None):
        if indices is None:
            indices=self.regions[name].indices
        return self._anuga.domain.get_quantity("stage").get_values(location=location,indices=indices)
    def get_elevation_ANUGA(self, name=None):
        return self.elevation[name]
    def get_rate_operator(self, name=None):
        return self.rate_operator[name]
    def get_area_ANUGA(self, name=None,indices=None,location="centroids"):
        if indices is None:
            indices=self.regions[name].indices
        return self._anuga.domain.areas[indices]



