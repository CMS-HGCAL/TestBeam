#--------------------------------------------------------------------------
# Basic assumption: a test beam detector can be modeled as a sequence of
# elements along the z axis (beam axis).
#--------------------------------------------------------------------------
# Components
#
# This lists all components from which a test beam geometry can be built.
# (x, y, z) is the location of the center of the component, which will be
# automatically determined given the location of the first component, as
# specified in the Geometry block. For sensitive components, set the
# sensitive attribute to True. A component is modeled as a Python dictionary, 
# a construct comprising a key and a value. The key is the name of the 
# component and the value is itself a dictionary containing key, value pairs 
# that define the attributes of the component.
#
# Every component must have (at a minimum) the following attributes:
#
#   shape, material, units (of length), side, thickness, (x, y, z).
# 
# As noted, the center of the component (x, y, z) is determined automatically.
# The basic assumption is that a test beam detector is simply a sequence of
# components along the beam axis (i.e., the z axis).
#
#
# World
#
# This specifies the world volume
#
# Geomtry
#
# This specifies the geometry, which is built from components and 
# composite components. The latter is a (Python) list of components.
#
#
# WARNING: Be careful about the placement of commas! It is easy to leave one
# out and confuse Python.
#--------------------------------------------------------------------------
Components = {'W21': {'shape': 'square',   # 2.1mm of W
                      'material': 'W',
                      'units': 'mm',
                      'side': 140.0,
                      'thickness': 2.1,
                      'x': 0.0,
                      'y': 0.0,
                      'z': 0.0},
              
              'W42':  {'shape': 'square',
                       'material': 'W',
                       'units': 'mm',
                       'side': 140.0,
                       'thickness': 4.2,
                       'x': 0.0,
                       'y': 0.0,
                       'z': 0.0},

              'Cu60': {'shape': 'hexagon',
                       'material': 'Cu',
                       'units': 'mm',
                       'side': 71.4598,
                       'thickness': 6.0,
                       'x': 0.0,
                       'y': 0.0,
                       'z': 0.0},

              'WCu06': {'shape': 'hexagon',
                        'material': 'WCu',
                        'units': 'mm',
                        'side': 71.4598,
                        'thickness': 0.6,
                        'x': 0.0,
                        'y': 0.0,
                        'z': 0.0},

              'Air30': {'shape': 'square',
                        'material': 'Air',
                        'units': 'mm',
                        'side': 71.4598,
                        'thickness': 3.0,
                        'x': 0.0,
                        'y': 0.0,
                        'z': 0.0},

              'Air60': {'shape': 'square',
                        'material': 'Air',
                        'units': 'mm',
                        'side': 71.4598,
                        'thickness': 6.0,
                        'x': 0.0,
                        'y': 0.0,
                        'z': 0.0},
              
              'Kapton':{'shape': 'hexagon',
                        'material': 'Air',
                        'units': 'mm',
                        'side': 71.4598,
                        'thickness': 0.01,
                        'x': 0.0,
                        'y': 0.0,
                        'z': 0.0},

              'Si020': {'shape': 'hexagon',
                        'material': 'Si',
                        'units': 'mm',
                        'sensitive': True,
                        'side': 71.4598,
                        'cellsize': 6.496345,
                        'thickness': 0.20,
                        'x': 0.0,
                        'y': 0.0,
                        'z': 0.0},

              'Si012': {'shape': 'hexagon',
                        'material': 'Si',
                        'units': 'mm',
                        'sensitive': False,
                        'side': 71.4598,
                        'thickness': 0.12,
                        'x': 0.0,
                        'y': 0.0,
                        'z': 0.0},

              # This is a composite component, built from the components
              # listed above. A composite component is modeled as an ordered
              # list of components
              'module2016_04': ['WCu06',
                                'Cu60',
                                'WCu06',
                                'Kapton',
                                'Si020',
                                'Si012'],

              'W42_Air60': ['W42',
                            'Air60'],

              'W21_Air60': ['W21',
                            'Air60'],

              'W21_Air30': ['W21',
                            'Air30']              
              }
#--------------------------------------------------------------------------
# World defines the overarching volume within which is placed a detector
# world volume, which in turn contains the detector.
#--------------------------------------------------------------------------
World = {'shape': 'box',
         'units': 'm',
         'xside': 0.5,
         'yside': 0.5,
         'zside': 1.0
         }

#--------------------------------------------------------------------------
# Geometry is specified as an ordered list of elements.
# An element is an instance of a component.
#--------------------------------------------------------------------------
# The Geometry block starts with a header containing the
# model, version number, and the desired location of the center of the
# first layer of the first element.
#
# The header block can also contain information specific to the model.

Geometry=[{'model':   5,
           'version': 1,
           'units': 'cm',
           'x':   0.0,   
           'y':   0.0,
           'z':  10.0
           },

          'W42_Air60',
          'W42_Air60',
          'W42_Air60',

          'W21_Air60',
          'W21_Air60',
          'W21_Air30',

          'module2016_04',
          'Air30',
          'module2016_04',
          'Air30',
          'module2016_04',
          'Air30',
          'module2016_04']
