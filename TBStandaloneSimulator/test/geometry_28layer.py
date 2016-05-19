from HGCal.TBStandaloneSimulator.Components import Components
#--------------------------------------------------------------------------
# Basic assumption: a test beam detector can be modeled as a sequence of
# elements along the z axis (beam axis).
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

Geometry=\
[{'model':   5,
  'version': 1,
  'units': 'cm',
  'x':   0.0,   
  'y':   0.0,
  'z':   5.0
  },
 
 'W2.0','Air1.0','module2.52','WCu0.6','Cu6.0','WCu0.6','module2.52','Air1.0',
 'W2.0','Air1.0','module2.52','WCu0.6','Cu6.0','WCu0.6','module2.52','Air1.0',
 'W2.0','Air1.0','module2.52','WCu0.6','Cu6.0','WCu0.6','module2.52','Air1.0',
 'W2.0','Air1.0','module2.52','WCu0.6','Cu6.0','WCu0.6','module2.52','Air1.0',
 'W2.0','Air1.0','module2.52','WCu0.6','Cu6.0','WCu0.6','module2.52','Air1.0',

 'W2.8','Air1.0','module2.52','WCu1.2','Cu6.0','WCu1.2','module2.52','Air1.0',
 'W2.8','Air1.0','module2.52','WCu1.2','Cu6.0','WCu1.2','module2.52','Air1.0',
 'W2.8','Air1.0','module2.52','WCu1.2','Cu6.0','WCu1.2','module2.52','Air1.0',
 'W2.8','Air1.0','module2.52','WCu1.2','Cu6.0','WCu1.2','module2.52','Air1.0',
 'W2.8','Air1.0','module2.52','WCu1.2','Cu6.0','WCu1.2','module2.52','Air1.0',

 'W4.2','Air1.0','module2.52','WCu2.2','Cu6.0','WCu2.2','module2.52','Air1.0',
 'W4.2','Air1.0','module2.52','WCu2.2','Cu6.0','WCu2.2','module2.52','Air1.0',
 'W4.2','Air1.0','module2.52','WCu2.2','Cu6.0','WCu2.2','module2.52','Air1.0',
 'W4.2','Air1.0','module2.52','WCu2.2','Cu6.0','WCu2.2','module2.52','Air1.0'
]

