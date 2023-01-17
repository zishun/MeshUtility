"""
This is part of MeshUtility

by Zishun Liu <liuzishun@gmail.com>
"""
'''
NOTE:
* MeshLab translation vector is special.
* FOV in Polyscope: ctrl+V does not change it. Set it in the menu bar.
'''


import numpy as np
import warnings
import xml.etree.ElementTree as ET
from scipy.spatial.transform import Rotation as Rotation


__all__ = ['parse_meshlab_view', 'generate_meshlab_view', 
           'parse_polyscope_view', 'generate_polyscope_view',
           'parse_opencv_matrices', 'generate_opencv_matrices']


def orthogonalize(R):
    '''
    Rotation matrix load from string should be orthogonalized.
    Otherwise, np.all(np.isclose(R.T @ R @ t, t)) is False.
    '''
    u, s, v = np.linalg.svd(R)
    return u@v


def focal_to_fovy(f, height):
    fovy = 2*np.arctan(height/2/f)*180/np.pi
    return fovy
    
    
def fovy_to_focal(fovy, height):
    f = height/2/np.tan(fovy/180*np.pi/2)
    return f


def parse_meshlab_view(xml_str):
    '''
    Reference: 
    https://github.com/cnr-isti-vclab/meshlab/blob/3c05b28ab83e0fd39ec732ca0aa629ef9475bdaf/src/meshlab/glarea.cpp
    https://github.com/cnr-isti-vclab/meshlab/blob/3c05b28ab83e0fd39ec732ca0aa629ef9475bdaf/src/meshlab/glarea.h
    https://github.com/cnr-isti-vclab/meshlab/blob/3c05b28ab83e0fd39ec732ca0aa629ef9475bdaf/src/common/utilities/load_save.cpp
    '''
    project = ET.fromstring(xml_str)
    for group in project:
        # print(group.tag)
        if group.tag == 'VCGCamera':
            R = np.fromstring(group.attrib['RotationMatrix'], dtype=float, sep=' ').reshape(4,4)
            canvas_size = np.fromstring(group.attrib['ViewportPx'], dtype=int, sep=' ')
            camera_type = int(group.attrib['CameraType'])
            t = np.fromstring(group.attrib['TranslationVector'], dtype=float, sep=' ')
            f = float(group.attrib['FocalMm'])
            pixel_size = np.fromstring(group.attrib['PixelSizeMm'], dtype=float, sep=' ')
            fov = focal_to_fovy(f/pixel_size[1], canvas_size[1])
            R[:3, :3] = orthogonalize(R[:3, :3])
            R[:3, 3] = R[:3,:3] @ t[:3]
        elif group.tag == 'ViewSettings':
            far = float(group.attrib['FarPlane'])
            near = float(group.attrib['NearPlane'])
        else:
            warnings.warn('Unknown group tag: %s'%group.tag)
    camera_types = ['perspective', 'orthographic']
    return {'type': camera_types[camera_type], 'width': canvas_size[0], 'height': canvas_size[1], 'fov': fov, 'near': near, 'far': far, 'view': R}
    
    # d = np.linalg.norm(t[:3])
    # y = np.tan(fov/2) * d
    # x = y / canvas_size[1] * canvas_size[0]


def generate_meshlab_view(cam):
    width = cam['width']
    height = cam['height']
    pixel_size = 0.036916077  # magic number in https://github.com/cnr-isti-vclab/meshlab/blob/3c05b28ab83e0fd39ec732ca0aa629ef9475bdaf/src/meshlab/glarea.cpp#L1997
    R = np.eye(4)
    R[:3,:3] = orthogonalize(cam['view'][:3,:3])
    R_str = np.array2string(R.ravel(), separator=' ')[1:-1]  # remove '[' and ']'
    t = np.ones((4,))
    t[:3] = R[:3,:3].T @ cam['view'][:3, 3]
    t_str = np.array2string(t.ravel(), separator=' ')[1:-1]  # remove '[' and ']'
    fmm = fovy_to_focal(cam['fov'], height)*pixel_size
    camera_type = 0 if cam['type']=='perspective' else 1
    result = '''<!DOCTYPE ViewState>
<project>
 <VCGCamera CameraType="%d" ViewportPx="%d %d" CenterPx="%d %d" BinaryData="0" LensDistortion="0 0" PixelSizeMm="%f %f" FocalMm="%f" RotationMatrix="%s" TranslationVector="%s"/>
 <ViewSettings NearPlane="%f" FarPlane="%f" TrackScale="1.0"/>
</project>
''' % (camera_type, width, height, width//2, height//2, pixel_size, pixel_size,
       fmm, R_str, t_str, cam['near'], cam['far']
      )
    return result
    
    d = np.linalg.norm(t[:3])
    fov = np.arctan(cam['top']/d) * 2
    fmm = fovy_to_focal(fov, height)*pixel_size


def parse_polyscope_view(dict_str, width=1280, height=720):
    '''
    Examples:
    '{"farClipRatio":20.0,"fov":45.0,"nearClipRatio":0.005,"projectionMode":"Perspective","viewMat":[0.350662589073181,-0.936452269554138,0.00963005423545837,-0.00613140314817429,0.336915373802185,0.135742008686066,0.931698441505432,-0.943150043487549,-0.873798370361328,-0.323467284440994,0.363104850053787,-4.57857847213745,0.0,0.0,0.0,1.0]}'
    '{"farClipRatio":20.0,"fov":45.0,"nearClipRatio":0.005,"projectionMode":"Orthographic","viewMat":[0.134992673993111,-0.990766167640686,0.0126598700881004,2.92058612103574e-05,0.293025732040405,0.0521232560276985,0.954683542251587,-0.281530499458313,-0.946528673171997,-0.125165998935699,0.297356098890305,-4.23833131790161,0.0,0.0,0.0,1.0]}'
    '''
    view_dict_polyscope = eval(dict_str)
    camera_type = view_dict_polyscope['projectionMode'].lower()
    fov = view_dict_polyscope['fov']
    near = view_dict_polyscope['nearClipRatio']
    far = view_dict_polyscope['farClipRatio']
    view = np.reshape(np.array(view_dict_polyscope['viewMat']), (4,4), order='C')
    return {'type': camera_type, 'fov': fov, 'width': width, 'height': height, 'near': near, 'far': far, 'view': view}
    
    d = np.linalg.norm(view[:3, 3])
    y = np.tan(fov/2) * d
    x = y / height * width


def generate_polyscope_view(cam):
    camera_type = cam['type'][0].upper()+cam['type'][1:]
    view_str = np.array2string(cam['view'].ravel(), separator=',')
    result = '{"farClipRatio":%f,"fov":%f,"nearClipRatio":%f,"projectionMode":"%s","viewMat":%s}' % (
        cam['far'], cam['fov'], cam['near'], camera_type, view_str
    )
    return result

    d = np.linalg.norm(cam['view'][:3, 3])
    fov = np.arctan(cam['top']/d) * 2
    

def parse_opencv_matrices(cameraMatrix, rvec, tvec):
    assert False, 'not implemented'


def generate_opencv_matrices(cam):
    '''
    Reference:
    Rodrigues formula in OpenCV
    https://docs.opencv.org/4.x/d9/d0c/group__calib3d.html#ga61585db663d9da06b68e70cfbf6a1eac
    '''
    fov = cam['fov']
    width = cam['width']
    height = cam['height']
    f = fovy_to_focal(fov, height)

    cameraMatrix = np.array([
        [f, 0, width/2-0.5],
        [0, f, height/2-0.5],
        [0, 0, 1]])
    Rt = np.array([
        [1, 0, 0],
        [0, -1, 0],
        [0, 0, -1]
    ]) @ cam['view'][:3, :]  # coordinate systems are different between OpenGL and OpenCV
    R = Rotation.from_matrix(Rt[:3, :3])
    rvec = R.as_rotvec()
    tvec = Rt[:, 3]
    return cameraMatrix, rvec, tvec


if __name__ == '__main__':
    pass
