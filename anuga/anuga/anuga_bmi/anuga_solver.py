#! /usr/bin/env python
import warnings

import numpy as np

import anuga

class AnugaSolver(object):

    def __init__(self, params):
                
        self._domain_type = str(params['domain_type'])
        self._shape = tuple(params['shape'])
        self._size = tuple(params['size'])
        self._boundary_filename = str(params['boundary_filename'])
        self._elevation_filename = str(params['elevation_filename'])
        self._output_filename = str(params['output_filename'])
        self._time_step = float(params['output_timestep'])
        self._bdry_tags = dict(params['boundary_tags'])
        self._bdry_conditions = dict(params['boundary_conditions'])
        self._stored_quantities = dict(params['stored_quantities'])
        self._max_triangle_area = float(params['maximum_triangle_area'])
        self._initial_flow_depth = float(params['initial_flow_depth'])
        
        self._interior_poly_filename = str(params['interior_polygon_filename'])
        self._interior_poly_triangle_area = float(params['interior_polygon_triangle_area'])
        self._interior_regions = None
        
        self._use_sed_operator = bool(params['toggle_sediment_transport'])
        self._inflow_concentration = float(params['inflow_sediment_concentration'])
        self._initial_concentration = float(params['initial_sediment_concentration'])
        self._use_veg_operator = bool(params['toggle_vegetation_drag'])
        self._veg_stem_diameter = params['vegetation_stem_diameter']
        self._veg_stem_spacing = params['vegetation_stem_spacing']
        self._mannings_n = float(params['Mannings_n_parameter'])
        
        self._elevation_profile = str(params['elevation_profile'])

        
        self._time = 0
        
        self.initialize_domain()
        self.set_boundary_conditions()
        self.set_other_domain_options()
        self.initialize_operators()
        
        
        # store initial elevations for differencing
        self._land_surface__initial_elevation = np.zeros_like(self.land_surface__elevation)
        
        
    def initialize_operators(self):
        """
        Initialize anuga operators
        
        Available operators:
        * Sed transport
        * Vegetation
        
        Operators remaining to implement:
        * Infiltration
        * Rainfall
        
        """
        
        if self._use_veg_operator:
            # better to have veg before sed transport??
            
            from anuga.operators.vegetation_operator import Vegetation_operator
            veg_op = Vegetation_operator(self.domain)
            
            self.land_vegetation__stem_spacing = self._veg_stem_spacing
            self.land_vegetation__stem_diameter = self._veg_stem_diameter
            
            
        
        if self._use_sed_operator:
            
            # this is the default flow algorithm but it's required for sed transport
            self.domain.set_flow_algorithm('DE0')
            
            from anuga.operators.sed_transport_operator import Sed_transport_operator
            sed_op = Sed_transport_operator(self.domain)
            
            
            assert self._initial_concentration <= 0.3, (
                    "Volumetric suspended sediment concentration must be <= 0.3")
                    
            assert self._inflow_concentration <= 0.3, (
                    "Inflow volumetric suspended sediment concentration must be <= 0.3")
            
            self.land_surface_water_sediment_suspended__volume_concentration = self._initial_concentration
            sed_op.set_inflow_concentration(self._inflow_concentration)
        
        
        
        
    def initialize_domain(self):
        """Initialize anuga domain"""


        assert self._domain_type[:5] in ['squar', 'recta', 'outli', 'irreg', 'bound'], (
            "Domain shape must be 'square'/'rectangle'/'rectangular' or "
            "'outline'/'irregular'. Shape '%s' is not recognized." %
            self._domain_type)



        if self._domain_type[:5] in ['squar', 'recta']:
        
            self.domain = anuga.rectangular_cross_domain(
                                self._shape[0],
                                self._shape[1],
                                len1 = self._size[0],
                                len2 = self._size[1])
                                
            # set some default values for quantities
            self.set_elevation_rectangular()
            
            
            
                                
        elif self._domain_type[:5] in ['outli', 'irreg', 'bound']:
        
            bounding_polygon = anuga.read_polygon(self._boundary_filename)
            
            filename_root = self._elevation_filename[:-4]
            
            assert self._elevation_filename[-4:] in ['.asc', '.pts'], (
                "Cannot recognize type of elevation file '%s'. "
                "Please use an .asc or .pts file." % self._elevation_filename)
            
            if self._elevation_filename[-4:] == '.asc':
                anuga.asc2dem(filename_root + '.asc')
                anuga.dem2pts(filename_root + '.dem')
            
            
            if ((self._interior_regions is None) and
                (self._interior_poly_triangle_area > 0.0) and
                (self._interior_poly_filename is not None)):
                
                try:
                    interior_poly = anuga.read_polygon(
                                    self._interior_poly_filename)
                    self._interior_regions = [[interior_poly,
                                             self._interior_poly_triangle_area]]
                except:
                    pass
            
                               
            # generalize the mesh creation to be able to use the sed transport operator
            evolved_quantities =  ['stage', 'xmomentum', 'ymomentum', 'concentration']
            
            anuga.pmesh.mesh_interface.create_mesh_from_regions(
                                bounding_polygon = bounding_polygon,
                                boundary_tags = self._bdry_tags,
                                interior_regions = self._interior_regions,
                                maximum_triangle_area = self._max_triangle_area,
                                filename = filename_root + '.msh')
                                
            self.domain = anuga.Domain(filename_root + '.msh',
                                       evolved_quantities = evolved_quantities)
            
            
            self.domain.set_quantity('elevation',
                                     filename = filename_root + '.pts')
            
        
        self.land_surface_water__depth = self._initial_flow_depth
        self.land_surface_water_surface__elevation += self._initial_flow_depth
        self.manning_n_parameter = self._mannings_n
                                 
        anuga.Quantity(self.domain, name='veg_diameter', register=True)
        anuga.Quantity(self.domain, name='veg_spacing', register=True)
        anuga.Quantity(self.domain, name='shear_stress', register=True)
        anuga.Quantity(self.domain, name='concentration', register=True)
        


    def set_elevation_rectangular(self):
    
            if self._elevation_profile == 'shallow linear ramp':
            
                self.domain.set_quantity('elevation', lambda x,y: -x/5.)
            
            if self._elevation_profile == 'steep linear ramp':
            
                self.domain.set_quantity('elevation', lambda x,y: -x/10.)

            if self._elevation_profile == 'tall step down':
            
                z = -x/10
                N = len(x)
                minx = np.floor(np.max(x)/2)
                maxx = np.ceil(np.max(x)/2 + 1)
                stepx = np.min(x[(x >= minx) & (x <= maxx)])
                stepz = 0.6 * (np.max(z) - np.min(z))

                for i in range(N):
                    if x[i] > stepx:
                        z[i] -= stepz
                        
                self.domain.set_quantity('elevation', z)            
            
            if self._elevation_profile == 'short step down':
            
                z = -x/10
                N = len(x)
                minx = np.floor(np.max(x)/2)
                maxx = np.ceil(np.max(x)/2 + 1)
                stepx = np.min(x[(x >= minx) & (x <= maxx)])
                stepz = 0.2 * (np.max(z) - np.min(z))

                for i in range(N):
                    if x[i] > stepx:
                        z[i] -= stepz
                        
                self.domain.set_quantity('elevation', z)
            
            
            if self._elevation_profile == 'step up':
            
                z = -x/10
                N = len(x)
                minx = np.floor(np.max(x)/2)
                maxx = np.ceil(np.max(x)/2 + 1)
                stepx = np.min(x[(x >= minx) & (x <= maxx)])
                stepz = 0.3 * (np.max(z) - np.min(z))

                for i in range(N):
                    if x[i] > stepx:
                        z[i] += stepz
                        
                self.domain.set_quantity('elevation', z)
            
            
            if self._elevation_profile == 'dam':
            
                z = -x/10.
                N = len(x)
                minx = np.floor(np.max(x)/2)
                maxx = np.ceil(np.max(x)/2)
                maxz = 0.9 * np.max(z)
    
                if len(x[(x >= minx) & (x <= maxx)]) < 5:
                    maxx = np.ceil(np.max(x)/2 + 1)
    
                for i in range(N):
                    if minx <= x[i] <= maxx:
                        z[i] = maxz
                        
                self.domain.set_quantity('elevation', z)
                       
                       
            if self._elevation_profile == 'thin cylinder':
                                              
                z = -x/10
                N = len(x)

                minx = np.floor(np.max(x)/2)
                maxx = np.ceil(np.max(x)/2 + 1)
                stepx = np.min(x[(x >= minx) & (x <= maxx)])

                miny = np.floor(np.max(y)/2)
                maxy = np.ceil(np.max(y)/2 + 1)
                stepy = np.min(y[(y >= miny) & (y <= maxy)])

                radius = 0.1 * (np.max(y) - np.min(y))

                for i in range(N):
                    if (x[i]-stepx)**2 + (y[i]-stepy)**2 < radius**2:
                        z[i] = np.max(z)
                        
                self.domain.set_quantity('elevation', z)
            
            
            if self._elevation_profile == 'thick cylinder':
                                              
                z = -x/10
                N = len(x)

                minx = np.floor(np.max(x)/2)
                maxx = np.ceil(np.max(x)/2 + 1)
                stepx = np.min(x[(x >= minx) & (x <= maxx)])

                miny = np.floor(np.max(y)/2)
                maxy = np.ceil(np.max(y)/2 + 1)
                stepy = np.min(y[(y >= miny) & (y <= maxy)])

                radius = 0.4 * (np.max(y) - np.min(y))

                for i in range(N):
                    if (x[i]-stepx)**2 + (y[i]-stepy)**2 < radius**2:
                        z[i] = np.max(z)
                        
                self.domain.set_quantity('elevation', z)
                
                
                
            if self._elevation_profile == 'wing dams':
            
                z = -x/10
                N = len(x)
                minx = np.floor(np.max(x)/2)
                stepx1 = np.min(x[(x >= minx)])
                stepx2 = np.min(x[(x > stepx1)])
                dist = 0.2 * (np.max(y) - np.min(y))
    
                for i in range(N):
                    if stepx1 <= x[i] <= stepx2:
                        if (y[i] < dist) or (y[i] > np.max(y) - dist):
                            z[i] += 1
                            
                self.domain.set_quantity('elevation', z)
                
                
            if self._elevation_profile == 'alternating dykes':
            
                z = -x/10
                N = len(x)

                minx = np.floor(np.max(x)/4)
                stepx1 = np.min(x[(x >= minx)])
                stepx2 = np.min(x[(x > stepx1)])

                minx = np.floor(np.max(x)/2)
                stepx3 = np.min(x[(x >= minx)])
                stepx4 = np.min(x[(x > stepx3)])

                minx = np.floor(3*np.max(x)/4)
                stepx5 = np.min(x[(x >= minx)])
                stepx6 = np.min(x[(x > stepx5)])

                dist = 0.3 * (np.max(y) - np.min(y))

                for i in range(N):
                    if stepx1 <= x[i] <= stepx2:
                        if (y[i] < dist):
                            z[i] += 1
            
                    if stepx3 <= x[i] <= stepx4:
                        if (y[i] > np.max(y) - dist):
                            z[i] += 1
            
                    if stepx5 <= x[i] <= stepx6:
                        if (y[i] < dist):
                            z[i] += 1
                            
                self.domain.set_quantity('elevation', z)
        
        
        
    def set_boundary_conditions(self):
        """
        Set boundary conditions for ANUGA domain
        
        Valid boundaries are (case insensitive):
        - reflective
        - transmissive
        - dirichlet / fixed (must specify stage at this boundary)
        - time (need to specify a lambda function)
        
        TODO:
        - check possible failure modes (how would anuga normally fail if the boundaries
            are not specified correctly?)
        - add defaults if not enough boundaries are specified
        - fail if receives an unknown boundary type
        - accept other inputs for time boundary (file?)
        """
        
        _bdry_conditions = {}
        
        for key, value in self._bdry_conditions.items():

            if isinstance(value, str):
                value = [value]
        
            bdry_type = value[0]
            
            if bdry_type.lower() == 'reflective':
                _bdry_conditions[key] = anuga.Reflective_boundary(self.domain)
                
                
                
            elif bdry_type.lower() == 'transmissive':
                _bdry_conditions[key] = anuga.Transmissive_boundary(self.domain)
                
                
                
            elif bdry_type.lower() in ['dirichlet', 'fixed']:
            
                dirichlet_vals = value[1:]
                
                assert len(value) > 1, ("Need to specify stage of "
                                        "Dirichlet boundary '%s'" % key)
                
                if len(dirichlet_vals) == 1: dirichlet_vals += [0., 0.]
                    
                _bdry_conditions[key] = anuga.Dirichlet_boundary(dirichlet_vals)
                
                
                
            elif bdry_type.lower() == 'time':
            
                assert len(value) > 1, ("Need to specify lambda function for "
                                        "Time boundary '%s'" % key)
                
                _bdry_conditions[key] = anuga.Time_boundary(domain = domain,
                                                                 function = value[1])
                
            else:
            
                raise ValueError("Did not recognize boundary type '%s' "
                               "of boundary '%s'" % (bdry_type, key))
                
        
        self._bdry_conditions = _bdry_conditions        
        self.domain.set_boundary(self._bdry_conditions)
        

    def set_other_domain_options(self):
    
        self.domain.set_name(self._output_filename)
        self.domain.set_quantities_to_be_stored(self._stored_quantities)                
        

    @property
    def grid_x(self):
        """x position of centroids"""
        var_values = self.domain.quantities['x'].centroid_values
        return var_values
        
    @property
    def grid_y(self):
        """y position of centroids"""
        var_values = self.domain.quantities['y'].centroid_values
        return var_values
        
    @property
    def grid_z(self):
        """z position of centroids"""
        return self.land_surface__elevation()

    @property
    def time_step(self):
        """The time step."""
        return self._time_step

    @time_step.setter
    def time_step(self, new_dt):
        self._time_step = new_dt

    @property
    def shape(self):
        """Number of grid rows and columns."""
        return self._shape

    @property
    def size(self):
        """Size of grid in meters."""
        return self._size


    ########
    
    @property
    def manning_n_parameter(self):
        """Manning's friction parameter"""
        return self.domain.quantities['friction'].centroid_values
        
    @manning_n_parameter.setter
    def manning_n_parameter(self, new_friction):
        self.domain.set_quantity('friction', new_friction, location='centroids')

    @property
    def land_surface__elevation(self):
        return self.domain.quantities['elevation'].centroid_values
        
    @land_surface__elevation.setter
    def land_surface__elevation(self, new_elev):
        self.domain.set_quantity('elevation', new_elev, location='centroids')

    @property
    def land_surface_water_surface__elevation(self):
        """Temperature values on the grid."""
        return self.domain.quantities['stage'].centroid_values

    @land_surface_water_surface__elevation.setter
    def land_surface_water_surface__elevation(self, new_stage):
    
        self.domain.set_quantity('stage', new_stage, location='centroids')
        
        new_depth = new_stage - self.land_surface__elevation
        self.domain.set_quantity('height', new_depth, location='centroids')
        
        
    @property
    def land_surface_water__depth(self):
        return self.domain.quantities['height'].centroid_values
        
    @land_surface_water__depth.setter
    def land_surface_water__depth(self, new_depth):
    
        self.domain.set_quantity('height', new_depth, location='centroids')
        
        new_stage = self.land_surface__elevation + new_depth
        self.domain.set_quantity('stage', new_stage, location='centroids')
        
        
    @property
    def land_surface_water_flow__x_component_of_momentum(self):
        return self.domain.quantities['xmomentum'].centroid_values
        
    @land_surface_water_flow__x_component_of_momentum.setter
    def land_surface_water_flow__x_component_of_momentum(self, new_mom):
        self.domain.set_quantity('xmomentum', new_mom, location='centroids')
        
    @property
    def land_surface_water_flow__y_component_of_momentum(self):
        return self.domain.quantities['ymomentum'].centroid_values
        
    @land_surface_water_flow__y_component_of_momentum.setter
    def land_surface_water_flow__y_component_of_momentum(self, new_mom):
        self.domain.set_quantity('ymomentum', new_mom, location='centroids')
        
    @property
    def land_surface_water_flow__shear_stress(self):
        return self.domain.quantities['shear_stress'].centroid_values
        
    @land_surface_water_flow__shear_stress.setter
    def land_surface_water_flow__shear_stress(self, new_ss):
        self.domain.set_quantity('shear_stress', new_ss, location='centroids')
    
    @property
    def land_surface_water_sediment_suspended__volume_concentration(self):
        return self.domain.quantities['concentration'].centroid_values
        
    @land_surface_water_sediment_suspended__volume_concentration.setter
    def land_surface_water_sediment_suspended__volume_concentration(self, new_c):
        self.domain.set_quantity('concentration', new_c, location='centroids')
        
    @property
    def land_vegetation__stem_spacing(self):
        return self.domain.quantities['veg_spacing'].centroid_values
        
    @land_vegetation__stem_spacing.setter
    def land_vegetation__stem_spacing(self, new_vs):
    
        if (type(new_vs) == str):
            self.domain.set_quantity('veg_spacing', filename = new_vs)
        else:
            self.domain.set_quantity('veg_spacing', new_vs, location='centroids')
        
    @property
    def land_vegetation__stem_diameter(self):
        return self.domain.quantities['veg_diameter'].centroid_values
        
    @land_vegetation__stem_diameter.setter
    def land_vegetation__stem_diameter(self, new_vd):
    
        if (type(new_vd) == str):
            self.domain.set_quantity('veg_diameter', filename = new_vd)
        else:
            self.domain.set_quantity('veg_diameter', new_vd, location='centroids')

    #########
    
    @property
    def land_surface__initial_elevation(self):
        """Initial surface elevation, for differencing"""
        return self._land_surface__initial_elevation
        
    @land_surface__initial_elevation.setter
    def land_surface__initial_elevation(self, new_init_elev):
        """Able to set the initial land surface elevation to
        set an initial elev for differencing. Use:
            
        land_surface__initial_elevation = land_surface__elevation
        """
        self._land_surface__initial_elevation[:] = new_init_elev
        
        
    def update_elev_difference(self):
        return self.land_surface__elevation - self._land_surface__initial_elevation


    #########    

    def update(self):
        """Evolve."""
        
        for t in self.domain.evolve(yieldstep = self._time_step, finaltime = self._time):
            print(self.domain.timestepping_statistics())

