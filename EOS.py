import woma
import matplotlib.pyplot as plt
import numpy as np
from contourpy import contour_generator
from scipy.interpolate import CubicSpline

path_to_SESAME_table = 'AQUA_H20.txt'

woma.load_eos_tables(['AQUA', 'CD21_HHe'])

def make_adiabat(P, T, mat_id):

    if mat_id == 304:
        A1_rho = woma.sesame.A1_rho_AQUA
        A1_T = woma.sesame.A1_T_AQUA
        A2_s = woma.sesame.A2_s_AQUA

    elif mat_id == 307:
        A1_rho = woma.sesame.A1_rho_CD21_HHe
        A1_T = woma.sesame.A1_T_CD21_HHe
        A2_s = woma.sesame.A2_s_CD21_HHe

    s_contour_gen = contour_generator(A1_T, A1_rho, A2_s)

    rho = woma.rho_P_T(P, T, 304)
    s = woma.s_rho_T(rho, T, 304)

    print(f'Making entropy contour at P = {P:.2e}, T = {T:.0f}, s = {s:.0f}')

    s_contours = s_contour_gen.lines(s)
    T_range_max = 0

    for s_contour in s_contours:

        min_T, max_T = np.min(s_contour[:, 0]), np.max(s_contour[:, 0])
        T_range = max_T - min_T

        if T_range > T_range_max:

            s_contour_rho = s_contour[:, 1]
            s_contour_T = s_contour[:, 0]

            T_range_max = T_range

    s_contour_P = woma.A1_P_rho_T(s_contour_rho, s_contour_T, np.full_like(s_contour_rho, 304))
    P_mask = (1e-6 < s_contour_P) & (s_contour_P < 2e9)

    s_contour_P = s_contour_P[P_mask][::-1]
    s_contour_rho = s_contour_rho[P_mask][::-1]
    s_contour_T = s_contour_T[P_mask][::-1]

    plt.plot(s_contour_P, s_contour_rho)
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig('fig.png', dpi=1000)

    interpolator = CubicSpline(s_contour_P, s_contour_rho)

    return interpolator
