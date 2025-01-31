import numpy as np
from EOS import make_adiabat
from analyse import get_snapshot_data

class atmosphere:

    def __init__(self, file_path, n, R_surface):
        
        pos, m, id, rho, T, P, S, u, mat_id, vel = get_snapshot_data(n, file_path)

        x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]
        r = np.sqrt(x ** 2 + y ** 2 + z ** 2)

        surface_mask = (r < R_surface) & (r > R_surface * 0.9)

        np.sum()

