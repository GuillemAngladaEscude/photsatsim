import numpy as np
import aux_functions as auxf


class Optic:

    def __init__(self, N_pixels: int, fov: float):

        self.N_pix = N_pixels
        self.fov = fov*np.pi/180


    def check_if_inframe(self,r_SBFOVF):
        if np.absolute(r_SBFOVF[2]) < np.sin(self.fov) and np.absolute(r_SBFOVF[0]) < np.sin(self.fov):
            return True
        else:
            return False


    def project_SBFOVF_2_sensor(self, r_SBFOVF, method):
        if method=='linear':
            fov = self.fov
            N_pixels = self.N_pix

            r_SS_SBFOVF = [0,1,0]
            TL_corner = auxf.R1(fov/2) @ auxf.R3(fov/2) @ r_SS_SBFOVF
            TR_corner = auxf.R1(fov/2) @ auxf.R3(-fov/2) @ r_SS_SBFOVF
            BL_corner = auxf.R1(-fov/2) @ auxf.R3(fov/2) @ r_SS_SBFOVF
            BR_corner = auxf.R1(-fov/2) @ auxf.R3(-fov/2) @ r_SS_SBFOVF

            #print(f"TL_corner: {TL_corner}")
            #print(f"TR_corner: {TR_corner}")
            #print(f"BL_corner: {BL_corner}")
            #print(f"BR_corner: {BR_corner}")

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

