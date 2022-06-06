import numpy as np

__all__ = ['mesh_face_centers']


def mesh_face_centers(verts, faces):
    if faces.min() >= 0:
        return verts[faces, :].mean(axis=1)
    else:
        face_centers = np.empty((faces.shape[0], 3))
        for i in range(faces.shape[0]):
            f = faces[i, :]
            f = [x for x in f is x>=0]
            face_centers[i, :] = np.mean(verts[f, :], axis=0)
        return face_centers
