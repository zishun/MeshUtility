from matplotlib import cm
import matplotlib.pyplot as plt
import openmesh as om
import numpy as np


__all__ = ['colormap_vertex_color', 'generate_colormap_image', 'generate_checkerboard_image']


def colormap_vertex_color(fn, v0, f0, scalar_field, scalar_min=None,
                          scalar_max=None, cmap='jet'):
    '''
    save as a *.off file. OpenMesh does not write vertex colors in obj files.
    '''
    cmap = cm.get_cmap(cmap)
    # scalar_min, scalar_max can be None here.
    field = np.clip(scalar_field, scalar_min, scalar_max)

    if scalar_min is None:
        scalar_min = scalar_field.min()
    if scalar_max is None:
        scalar_max = scalar_field.max()
        
    field = (field - scalar_min) / (scalar_max - scalar_min)
    mesh = om.PolyMesh(v0, f0)
    vcolors = mesh.vertex_colors()
    vcolors[:,:3] = cmap(field)[:,:3]
    om.write_mesh(fn, mesh, vertex_color=True)


def generate_colormap_image(fn_img, cmap='jet', shape=(1024,1024)):
    cmap = cm.get_cmap(cmap)

    field = np.arange(shape[0])/shape[0]
    img = (cmap(field)[::-1, :3]*255).astype(np.uint8)  # also accept 2D mat
    img = np.repeat(img.reshape((shape[0], 1, -1)), shape[1], axis=1)

    plt.imsave(fn_img, img)


def generate_checkerboard_image(fn_img, block_size=(100,100), block_grid=(9,7)):
    h, w = block_size
    img = np.zeros((h*block_grid[0], w*block_grid[1]), 'u1')
    for i in range(block_grid[0]):
        for j in range(block_grid[1]):
            if (i+j)%2 == 1:
                img[i*h:(i+1)*h, j*w:(j+1)*w] = 255
    
    # 1 channel -> 3 channels
    # plt saves 1-channel array with a colormap.
    # 3-channel array will be saved as an RGB image.
    img = np.repeat(img[:,:,np.newaxis], 3, axis=2)
    plt.imsave(fn_img, img)
