import numpy as np
import swiftsimio as sw
from swiftsimio.visualisation.slice import slice_gas
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm, Normalize
import unyt
import h5py
import woma

R_earth = 6.371e6  # m
M_earth = 5.9724e24  # kg


woma.load_eos_tables(['ANEOS_forsterite', 'ANEOS_iron', 'AQUA'])

path = '/home/pt426/Impacts/snapshots'
# path = '/home/pavan/Impacts/sub-neptune-impact-2/snapshots'
prefix = 'sub_neptune_impact'

def plot(n, plot_T=False, plot_P=False, file_path=''):
    
    if file_path != '':
        data = sw.load(f'{file_path}_{n:04.0f}.hdf5')
    else:
        data = sw.load(f'{path}/{prefix}_{n:04.0f}.hdf5')

    data.gas.internal_energies.convert_to_mks()
    data.gas.densities.convert_to_mks()

    pos, m = np.array(data.gas.coordinates), np.array(data.gas.masses)

    COM = np.sum(pos.T * m, axis=1) / np.sum(m)

    u, rho, mat_id = np.array(data.gas.internal_energies), np.array(data.gas.densities), np.array(data.gas.material_ids)
    mat_id[mat_id == 901] = 304

    T = woma.A1_T_u_rho(u, rho, mat_id)
    P = woma.A1_P_u_rho(u, rho, mat_id)

    data.gas.temperatures = sw.objects.cosmo_array(T * unyt.K)
    data.gas.temperatures.cosmo_factor = data.gas.internal_energies.cosmo_factor

    data.gas.pressures = sw.objects.cosmo_array(P * unyt.Pa)
    data.gas.pressures.cosmo_factor = data.gas.internal_energies.cosmo_factor

    center = [COM[0], COM[1]]
    size = 3
    region = [center[0] - size, center[0] + size, center[1] - size, center[1] + size] * unyt.Rearth

    mass_map = slice_gas(
        data,
        z_slice=0.5 * data.metadata.boxsize[2],
        resolution=1024,
        project='masses',
        region=region,
        parallel=True
    )

    data.gas.mass_weighted_temperatures = data.gas.masses * data.gas.temperatures
    data.gas.mass_weighted_pressures = data.gas.masses * data.gas.pressures

    mass_weighted_temperatures_map = slice_gas(
        data,
        z_slice=0.5 * data.metadata.boxsize[2],
        resolution=1024,
        project="mass_weighted_temperatures",
        region=region,
        parallel=True,
    )

    mass_weighted_pressures_map = slice_gas(
        data,
        z_slice=0.5 * data.metadata.boxsize[2],
        resolution=1024,
        project="mass_weighted_pressures",
        region=region,
        parallel=True,
    )

    temperatures_map = mass_weighted_temperatures_map / mass_map
    pressures_map = mass_weighted_pressures_map / mass_map

    mass_map.convert_to_units(unyt.g / unyt.cm ** 3)
    temperatures_map.convert_to_units(unyt.K)
    pressures_map.convert_to_units(unyt.Pa)

    plt.imshow(mass_map.value, cmap='magma')
    plt.colorbar()
    plt.show()

    if plot_T:
        plt.imshow(temperatures_map.value, cmap='turbo', norm=Normalize(vmin=0, vmax=4000))
        plt.colorbar()
        plt.show()

    if plot_P:
        plt.imshow(pressures_map.value, cmap='plasma', norm=LogNorm())
        plt.colorbar()
        plt.show()

def get_snapshot_data(n, file_path=''):

    if file_path != '':
        data = sw.load(f'{file_path}_{n:04.0f}.hdf5')
    else:
        data = sw.load(f'{path}/{prefix}_{n:04.0f}.hdf5')

    data.gas.internal_energies.convert_to_mks()
    data.gas.densities.convert_to_mks()
    data.gas.masses.convert_to_mks()

    pos, m = np.array(data.gas.coordinates), np.array(data.gas.masses)
    vel = np.array(data.gas.velocities)
    id = np.array(data.gas.particle_ids)
    u, rho, mat_id = np.array(data.gas.internal_energies), np.array(data.gas.densities), np.array(data.gas.material_ids)

    COM = np.sum(pos.T * m, axis=1) / np.sum(m)
    pos -= COM

    mat_id[mat_id == 901] = 304

    T = woma.A1_T_u_rho(u, rho, mat_id)
    P = woma.A1_P_u_rho(u, rho, mat_id)

    return pos, m, id, rho, T, P, u, mat_id, vel
