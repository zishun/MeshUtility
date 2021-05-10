USE_CPP = True  # set False if those CPP-based modules are not needed

from .src.colormap import *
from .src.isocurve import *
from .src.mesh_cut import *
from .src.mesh_split import *
from .src.mesh2graph import *
from .src.sphere_cvt import *
from .src.write_obj_lines import *
if USE_CPP:
    from .src.lib import pygeodesic
    from .src.lib import pyisocurve
