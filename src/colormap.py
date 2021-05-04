from matplotlib import cm
import openmesh as om
import numpy as np

__all__ = ['colormap_vertex_color']


def colormap_vertex_color(fn, v0, f0, scalar_field, scalar_min=None,
                          scalar_max=None):
    if scalar_min is None and scalar_max is None:
        field = scalar_field.copy()
    else:
        field = np.clip(scalar_field, scalar_min, scalar_max)    
    m0 = field.min()
    m1 = field.max()
    field = (field - m0) / (m1 - m0)
    mesh = om.PolyMesh(v0, f0)
    # mesh.request_vertex_colors()
    for vh in mesh.vertices():
        c = cm.jet(field[vh.idx()])
        mesh.set_color(vh, c)
    om.write_mesh(fn, mesh, vertex_color=True)
