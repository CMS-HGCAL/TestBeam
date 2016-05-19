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
# Geometry
#
# This specifies the geometry, which is built from components and 
# composite components. The latter is a (Python) list of components.
#
#
# WARNING: Be careful about the placement of commas! It is easy to leave one
# out and confuse Python.
#--------------------------------------------------------------------------
Components = {'W2.0': {'shape': 'hexagon',
                       'material': 'W',
                       'units': 'mm',
                       'side': 71.4598,
                       'thickness': 2.0,
                       'x': 0.0,
                       'y': 0.0,
                       'z': 0.0},

              'W2.1': {'shape': 'hexagon',
                       'material': 'W',
                       'units': 'mm',
                       'side': 71.4598,
                       'thickness': 2.1,
                       'x': 0.0,
                       'y': 0.0,
                       'z': 0.0},

              'W2.8': {'shape': 'hexagon',
                       'material': 'W',
                       'units': 'mm',
                       'side': 71.4598,
                       'thickness': 2.8,
                       'x': 0.0,
                       'y': 0.0,
                       'z': 0.0},
              
              'W4.2':  {'shape': 'hexagon',
                        'material': 'W',
                        'units': 'mm',
                        'side': 71.4598,
                        'thickness': 4.2,
                        'x': 0.0,
                        'y': 0.0,
                        'z': 0.0},

              'Cu6.0': {'shape': 'hexagon',
                        'material': 'Cu',
                        'units': 'mm',
                        'side': 71.4598,
                        'thickness': 6.0,
                        'x': 0.0,
                        'y': 0.0,
                        'z': 0.0},

              'WCu0.6': {'shape': 'hexagon',
                         'material': 'WCu',
                         'units': 'mm',
                         'side': 71.4598,
                         'thickness': 0.6,
                         'x': 0.0,
                         'y': 0.0,
                         'z': 0.0},

              'WCu1.2': {'shape': 'hexagon',
                         'material': 'WCu',
                         'units': 'mm',
                         'side': 71.4598,
                         'thickness': 1.2,
                         'x': 0.0,
                         'y': 0.0,
                         'z': 0.0},

              'WCu2.2': {'shape': 'hexagon',
                        'material': 'WCu',
                        'units': 'mm',
                        'side': 71.4598,
                        'thickness': 2.2,
                        'x': 0.0,
                        'y': 0.0,
                        'z': 0.0},

              'Air1.0': {'shape': 'hexagon',
                         'material': 'Air',
                         'units': 'mm',
                         'side': 71.4598,
                         'thickness': 1.0,
                         'x': 0.0,
                         'y': 0.0,
                         'z': 0.0},

              'Air0.2': {'shape': 'hexagon',
                         'material': 'Air',
                         'units': 'mm',
                         'side': 71.4598,
                         'thickness': 0.2,
                         'x': 0.0,
                         'y': 0.0,
                         'z': 0.0},

              'Air3.0': {'shape': 'hexagon',
                         'material': 'Air',
                         'units': 'mm',
                         'side': 71.4598,
                         'thickness': 3.0,
                         'x': 0.0,
                         'y': 0.0,
                         'z': 0.0},

              'Air6.0': {'shape': 'hexagon',
                         'material': 'Air',
                         'units': 'mm',
                         'side': 71.4598,
                         'thickness': 6.0,
                         'x': 0.0,
                         'y': 0.0,
                         'z': 0.0},
              
              'Kapton0.01':{'shape': 'hexagon',
                            'material': 'Air',
                            'units': 'mm',
                            'side': 71.4598,
                            'thickness': 0.01,
                            'x': 0.0,
                            'y': 0.0,
                            'z': 0.0},

              'G101.2':{'shape': 'hexagon',
                        'material': 'G10',
                        'units': 'mm',
                        'side': 71.4598,
                        'thickness': 1.2,
                        'x': 0.0,
                        'y': 0.0,
                        'z': 0.0},

              'FR40.5': {'shape': 'hexagon',
                         'material': 'FR4',
                         'units': 'mm',
                         'side': 71.4598,
                         'thickness': 0.5,
                         'x': 0.0,
                         'y': 0.0,
                         'z': 0.0},

              'Epoxy0.5': {'shape': 'hexagon',
                           'material': 'Epoxy',
                          'units': 'mm',
                           'side': 71.4598,
                           'thickness': 0.5,
                           'x': 0.0,
                           'y': 0.0,
                           'z': 0.0},

              'Si0.2': {'shape': 'hexagon',
                        'material': 'Si',
                        'units': 'mm',
                        'sensitive': True,
                        'side': 71.4598,
                        'cellsize': 6.496345,
                        'thickness': 0.20,
                        'x': 0.0,
                        'y': 0.0,
                        'z': 0.0},

              'Si0.1': {'shape': 'hexagon',
                        'material': 'Si',
                        'units': 'mm',
                        'sensitive': False,
                        'side': 71.4598,
                        'thickness': 0.1,
                        'x': 0.0,
                        'y': 0.0,
                        'z': 0.0},

              # The following are composite components, built from the 
              # components listed above. A composite component is modeled 
              # as an ordered list of components

              'module2016_04': ['WCu0.6',
                                'Cu6.0',
                                'WCu0.6',
                                'Kapton0.01',
                                'Si0.2',
                                'Si0.1'],

              'module2.52': ['FR40.5',
                             'G101.2',
                             'Kapton0.01',
                             'Kapton0.01',
                             'Si0.2',
                             'Si0.1',
                             'Epoxy0.5'],
              
              'moduleV1.0': ['W4.2',
                             'Air6.0',
                             'W2.8',
                             'Air3.0',
                             'Cu6.0',
                             'WCu2.2',
                             'Air0.2',
                             'Si0.1',
                             'Si0.2'],

              'moduleV2.0': ['W4.2',
                             'Air6.0',
                             'W2.8',
                             'Air3.0',
                             'Cu6.0',
                             'WCu1.2',
                             'Air0.2',
                             'Si0.1',
                             'Si0.2']
              }
