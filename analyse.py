import numpy as np
import swiftsimio as sw
from swiftsimio import load
from swiftsimio.visualisation.slice import slice_gas
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm, Normalize
import unyt
import h5py
import woma

R_earth = 6.371e6  # m
M_earth = 5.9724e24  # kg


woma.load_eos_tables(['ANEOS_forsterite', 'ANEOS_iron', 'AQUA'])

path = '/home/pavan/Impacts/sub-neptune-impact/snapshots'

def plot(n):
    data = load(f'{path}/sub_neptune_impact_{n:04.0f}.hdf5')

    data.gas.internal_energies.convert_to_mks()
    data.gas.densities.convert_to_mks()

    u, rho, mat_id = np.array(data.gas.internal_energies), np.array(data.gas.densities), np.array(data.gas.material_ids)
    mat_id[mat_id == 901] = 304

    T = woma.A1_T_u_rho(u, rho, mat_id)
    P = woma.A1_P_u_rho(u, rho, mat_id)

    data.gas.temperatures = sw.objects.cosmo_array(T * unyt.K)
    data.gas.temperatures.cosmo_factor = data.gas.internal_energies.cosmo_factor

    data.gas.pressures = sw.objects.cosmo_array(P * unyt.Pa)
    data.gas.pressures.cosmo_factor = data.gas.internal_energies.cosmo_factor

    mass_map = slice_gas(
        data,
        z_slice=0.5 * data.metadata.boxsize[2],
        resolution=1024,
        project='masses',
        region=[47, 53, 47, 53] * unyt.Rearth,
        parallel=True
    )

    data.gas.mass_weighted_temperatures = data.gas.masses * data.gas.temperatures
    data.gas.mass_weighted_pressures = data.gas.masses * data.gas.pressures

    mass_weighted_temperatures_map = slice_gas(
        data,
        z_slice=0.5 * data.metadata.boxsize[2],
        resolution=1024,
        project="mass_weighted_temperatures",
        region=[47, 53, 47, 53] * unyt.Rearth,
        parallel=True,
    )

    mass_weighted_pressures_map = slice_gas(
        data,
        z_slice=0.5 * data.metadata.boxsize[2],
        resolution=1024,
        project="mass_weighted_pressures",
        region=[47, 53, 47, 53] * unyt.Rearth,
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

    plt.imshow(temperatures_map.value, cmap='turbo', norm=Normalize(vmin=0, vmax=2000))
    plt.colorbar()
    plt.show()

    plt.imshow(pressures_map.value, cmap='plasma', norm=LogNorm())
    plt.colorbar()
    plt.show()

# plot(0)

data = load(f'{path}/sub_neptune_impact_{n:04.0f}.hdf5')

data.gas.internal_energies.convert_to_mks()
data.gas.densities.convert_to_mks()

ids = np.arrau(data.gas.particle_ids)

u, rho, mat_id = np.array(data.gas.internal_energies), np.array(data.gas.densities), np.array(data.gas.material_ids)
mat_id[mat_id == 901] = 304

T = woma.A1_T_u_rho(u, rho, mat_id)
P = woma.A1_P_u_rho(u, rho, mat_id)

data.gas.temperatures = sw.objects.cosmo_array(T * unyt.K)
data.gas.temperatures.cosmo_factor = data.gas.internal_energies.cosmo_factor

data.gas.pressures = sw.objects.cosmo_array(P * unyt.Pa)
data.gas.pressures.cosmo_factor = data.gas.internal_energies.cosmo_factor

f1 = h5py.File(f'{path}/sub_neptune_impact_0000.hdf5')
f2 = h5py.File(f'{path}/sub_neptune_impact_0160.hdf5')

x_pos_initial = f1['GasParticles']['Coordinates'][()][:, 0]

impactor_ids = f1['GasParticles']['ParticleIDs'][()][x_pos_initial > 54]

f1.close()

print(f2['GasParticles'].keys())

pos = f2['GasParticles']['Coordinates'][()]
m = f2['GasParticles']['Masses'][()] / M_earth
mat_ids = f2['GasParticles']['MaterialIDs'][()]

rho = f2['GasParticles']['Densities'][()]
u = f2['GasParticles']['InternalEnergies'][()]

rho *= (R_earth ** 3 / M_earth)
u *= (R_earth ** -2)

plt.hist(np.log10(rho))
plt.show()

mat_ids[mat_ids == 901] = 304

P = woma.A1_P_u_rho(u, rho, mat_ids)
T = woma.A1_T_u_rho(u, rho, mat_ids)

COM = np.sum(pos.T * m, axis=1) / np.sum(m)
pos -= COM

x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]
r = np.sqrt(x ** 2 + y ** 2 + z ** 2)

final_ids = f2['GasParticles']['ParticleIDs'][()]

impactor_mask = np.isin(final_ids, impactor_ids)
z_mask = (z > -0.1) & (z < 0.1)
forsterite_mask = (mat_ids == 400)
iron_mask = (mat_ids == 401)
water_mask = (mat_ids == 304)

s = 1

plt.scatter(x[z_mask & water_mask & ~impactor_mask], y[z_mask & water_mask & ~impactor_mask], s=s, c='blue')
plt.scatter(x[z_mask & forsterite_mask & ~impactor_mask], y[z_mask & forsterite_mask & ~impactor_mask], s=s, c='grey')
plt.scatter(x[z_mask & iron_mask & ~impactor_mask], y[z_mask & iron_mask & ~impactor_mask], s=s, c='black')
plt.scatter(x[z_mask & forsterite_mask & impactor_mask], y[z_mask & forsterite_mask & impactor_mask], s=s, c='magenta')
plt.scatter(x[z_mask & iron_mask & impactor_mask], y[z_mask & iron_mask & impactor_mask], s=s, c='red')
plt.show()

plt.scatter(r[water_mask & z_mask], P[water_mask & z_mask])
plt.yscale('log')
plt.show()

plt.scatter(T[water_mask & z_mask], P[water_mask & z_mask], c=r[water_mask & z_mask])
plt.colorbar()
plt.yscale('log')
plt.show()