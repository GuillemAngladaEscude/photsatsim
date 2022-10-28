import numpy as np
import aux_functions as auxf


class Optic:
    """
        A class used to define the parameters of the optical system on-board and to calculate the position of
        any given vector pointing to the sky in the local frame projected to the sensor

        This clss uses the numpy module

        Attributes
        ----------
        mag_max: int
                maximum star magnitude detectable by the sensor. Although it could be a non-integer number, at the current
                version it needs to be in order to load the correct catalog file among the library of catalogs filtered
                up to an integer magnitude

        Methods
        -------
        check_if_inframe(r_SBFOVF)
            Check if the sky position defined by a given vector in the Satellite-Based, Field-Of-View-Fixed (SBFOVF) Frame
            lays inside the frame or not

        project_pqr_2_uv(r_pqr):
            Project a given 3D, unit vector in the pqr frame to the 2D uv frame, contained in the sensor plane and centered
            at its centre.

        uv_2_xy(u, v):
             Convert the position of a point in the sensor given in the uv frame to the xy frame
        """

    def __init__(self, N_pixels_hor: int, N_pixels_vert: int, d_hor: float, d_vert: float, focal_length: float, mag_max: int):
#N_pixels_hor, N_pixels_vert, d_hor, d_vert,
        """
        Initialise the class Optic

        :param N_pixels_hor: Number of horizontal pixels of the sensor
        :param N_pixels_vert: Number of vertical pixels of the sensor
        :param d_hor: Physical sensor width [cm]
        :param d_vert: Physical sensor height [cm]
        :param focal_length: Focal length of the equivalent optical system [cm]
        :param mag_max: maximum star magnitude detectable by the sensor []
"""
        self.__N_pix_hor = N_pixels_hor
        self.__N_pix_vert = N_pixels_vert
        self.__fov_hor = 2*np.arctan(d_hor/(2*focal_length))
        self.__fov_vert = 2*np.arctan(d_vert/(2*focal_length))
        self.__d_hor = d_hor
        self.__d_vert = d_vert
        self.__f = focal_length
        self.mag_max = mag_max

    @property
    def fov_hor(self):
        return self.__fov_hor

    @property
    def fov_vert(self):
        return self.__fov_vert

    def check_if_inframe(self,r_SBFOVF):
        """
        Check if the sky position defined by a given vector in the Satellite-Based, Field-Of-View-Fixed (SBFOVF) Frame
        lays inside the frame or not.
        :param r_SBFOVF: unit vector pointing to the desired location in the Satellite-Based, Field-Of-View-Fixed
            (SBFOVF) Frame [ , , ]
        :return: True if the point in the sky given by the input vector is inside the frame, and False otherwise
        """
        if np.absolute(r_SBFOVF[2]) < np.sin(self.__fov_vert/2) and np.absolute(r_SBFOVF[0]) < np.sin(self.__fov_hor/2):
            return True
        else:
            return False

    def project_pqr_2_uv(self, r_pqr):
        """
        Project a given 3D, unit vector in the pqr frame to the 2D uv frame, contained in the sensor plane and centered
            at its centre.
        :param r_pqr: Unit vector in the pqr frame [ , , ]
        :return: Position, in the uv frame, of the point where r_pqr points [rad, rad]
        """
        # todo: aquí es tindran en compte els efectes de distorsió de l'òptica
        u, v = r_pqr[0], r_pqr[1]
        return [u, v]

    def uv_2_xy(self, u, v):
        """
        Convert the position of a point in the sensor given in the uv frame to the xy frame
        :param u: Horizontal component of the sensor point in the uv frame [rad]
        :param v: Vertical component of the sensor point in the uv frame [rad]
        :return: Position if the sensor point in the xy frame [pix, pix]
        """
        x = u * self.__f * self.__N_pix_hor/self.__d_hor
        y = v * self.__f * self.__N_pix_vert/self.__d_vert
        return [x, y]


    #old:
    # def project_SBFOVF_2_sensor(self, r_SBFOVF, method):
    #     if method=='linear':
    #         fov = self.fov
    #         N_pixels = self.N_pix
    #
    #         r_SS_SBFOVF = [0,1,0]
    #         TL_corner = auxf.R1(fov/2) @ auxf.R3(fov/2) @ r_SS_SBFOVF
    #         TR_corner = auxf.R1(fov/2) @ auxf.R3(-fov/2) @ r_SS_SBFOVF
    #         BL_corner = auxf.R1(-fov/2) @ auxf.R3(fov/2) @ r_SS_SBFOVF
    #         BR_corner = auxf.R1(-fov/2) @ auxf.R3(-fov/2) @ r_SS_SBFOVF
    #
    #         zmin = BL_corner[2]
    #         zmax = TL_corner[2]
    #         z = r_SBFOVF[2]
    #
    #         x_frame = N_pixels*(1 - (z-zmin)/(zmax-zmin))
    #
    #         xmin = BL_corner[0]
    #         xmax = BR_corner[0]
    #         x = r_SBFOVF[0]
    #
    #         y_frame = N_pixels*(x-xmin)/(xmax-xmin)
    #
    #         return [x_frame, y_frame]
    #
    #
    #     else:
    #         print("Projection method not acknowledged")

