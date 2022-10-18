import numpy as np
import aux_functions as auxf


class Optic:

    def __init__(self, N_pixels: int, fov: float, focal_length: float, mag_max: int):

        self.N_pix = N_pixels
        self.fov = fov*np.pi/180
        self.size = 2*focal_length*np.tan(fov*(np.pi/180)/2)
        self.f = focal_length
        self.mag_max = mag_max


    def check_if_inframe(self,r_SBFOVF):
        if np.absolute(r_SBFOVF[2]) < np.sin(self.fov/2) and np.absolute(r_SBFOVF[0]) < np.sin(self.fov/2):
            return True
        else:
            return False

    def project_pqr_2_uv(self, r_pqr):
        #todo: aquí es tindran en compte els efectes de distorsió de l'òptica
        u, v = r_pqr[0], r_pqr[1]
        return [u, v]

    def uv_2_xy(self, u, v):
        x = u * self.f * self.N_pix/self.size
        y = v * self.f * self.N_pix/self.size
        return [x, y]

    #old:
    def project_SBFOVF_2_sensor(self, r_SBFOVF, method):
        if method=='linear':
            fov = self.fov
            N_pixels = self.N_pix

            r_SS_SBFOVF = [0,1,0]
            TL_corner = auxf.R1(fov/2) @ auxf.R3(fov/2) @ r_SS_SBFOVF
            TR_corner = auxf.R1(fov/2) @ auxf.R3(-fov/2) @ r_SS_SBFOVF
            BL_corner = auxf.R1(-fov/2) @ auxf.R3(fov/2) @ r_SS_SBFOVF
            BR_corner = auxf.R1(-fov/2) @ auxf.R3(-fov/2) @ r_SS_SBFOVF

            zmin = BL_corner[2]
            zmax = TL_corner[2]
            z = r_SBFOVF[2]

            x_frame = N_pixels*(1 - (z-zmin)/(zmax-zmin))

            xmin = BL_corner[0]
            xmax = BR_corner[0]
            x = r_SBFOVF[0]

            y_frame = N_pixels*(x-xmin)/(xmax-xmin)

            return [x_frame, y_frame]


        else:
            print("Projection method not acknowledged")

