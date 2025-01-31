import numpy as np
import matplotlib.pyplot as plt
import woma

woma.load_eos_tables(['ANEOS_forsterite', 'ANEOS_iron', 'AQUA'])

R_earth = 6.371e6  # m
M_earth = 5.9724e24  # kg

def plot_spherical_profiles(planet):
    fig, ax = plt.subplots(2, 2, figsize=(8, 8))

    ax[0, 0].plot(planet.A1_r / R_earth, planet.A1_rho)
    ax[0, 0].set_xlabel(r"Radius, $r$ $[R_\oplus]$")
    ax[0, 0].set_ylabel(r"Density, $\rho$ [kg m$^{-3}$]")
    ax[0, 0].set_yscale("log")
    ax[0, 0].set_xlim(0, None)

    ax[1, 0].plot(planet.A1_r / R_earth, planet.A1_m_enc / M_earth)
    ax[1, 0].set_xlabel(r"Radius, $r$ $[R_\oplus]$")
    ax[1, 0].set_ylabel(r"Enclosed Mass, $M_{<r}$ $[M_\oplus]$")
    ax[1, 0].set_xlim(0, None)
    ax[1, 0].set_ylim(0, None)

    ax[0, 1].plot(planet.A1_r / R_earth, planet.A1_P)
    ax[0, 1].set_xlabel(r"Radius, $r$ $[R_\oplus]$")
    ax[0, 1].set_ylabel(r"Pressure, $P$ [Pa]")
    ax[0, 1].set_yscale("log")
    ax[0, 1].set_xlim(0, None)

    ax[1, 1].plot(planet.A1_r / R_earth, planet.A1_T)
    ax[1, 1].set_xlabel(r"Radius, $r$ $[R_\oplus]$")
    ax[1, 1].set_ylabel(r"Temperature, $T$ [K]")
    ax[1, 1].set_xlim(0, None)
    ax[1, 1].set_ylim(0, None)

    plt.tight_layout()
    plt.show()


sub_neptune = woma.Planet(
    name='Sub-Neptune',
    A1_mat_layer=["ANEOS_iron", "ANEOS_forsterite", "AQUA"],
    A1_T_rho_type=["adiabatic", "adiabatic", "adiabatic"],
    P_s=1e5,
    T_s=300,
    A1_R_layer=[0.5 * R_earth, 0.9 * R_earth, 2.3 * R_earth],
)

sub_neptune.gen_prof_L3_find_M_given_R_R1_R2(M_max=10 * M_earth)

mars = woma.Planet(
    name='Mars',
    A1_mat_layer=["ANEOS_iron", "ANEOS_forsterite"],
    A1_T_rho_type=["adiabatic", "adiabatic"],
    P_s=1e5,
    T_s=300,
    A1_R_layer=[0.25 * R_earth, 0.5 * R_earth]
)

mars.gen_prof_L2_find_M_given_R_R1(M_max=1 * M_earth)

moon = woma.Planet(
    name='Moon',
    A1_mat_layer=["ANEOS_iron", "ANEOS_forsterite"],
    A1_T_rho_type=["adiabatic", "adiabatic"],
    P_s=1e5,
    T_s=300,
    A1_R_layer=[350e3, 1700e3]
)

moon.gen_prof_L2_find_M_given_R_R1(M_max=1 * M_earth)

impactor = mars

n_particles = 1e6

M_total = sub_neptune.M + impactor.M

target_particles = woma.ParticlePlanet(sub_neptune, int(n_particles * (sub_neptune.M / M_total)), verbosity=1)
impactor_particles = woma.ParticlePlanet(impactor, int(n_particles * (impactor.M / M_total)), verbosity=1)

A1_pos_i, A1_vel_i = woma.impact_pos_vel_b_v_c_t(
    b=0.3,
    v_c=1.1,
    t=1800,
    R_t=sub_neptune.R,
    R_i=impactor.R,
    M_t=sub_neptune.M,
    M_i=impactor.M,
    units_v_c="v_esc",
    units_b="b"
)

impactor_particles.A2_pos += A1_pos_i
impactor_particles.A2_vel += A1_vel_i

A1_COM = np.sum(target_particles.A2_pos.T * target_particles.A1_m, axis=1) / np.sum(target_particles.A1_m)

pos = np.array(target_particles.A2_pos) - A1_COM

A1_vel_COM = (sub_neptune.M * A1_vel_i) / M_total

target_particles.A2_vel -= A1_vel_COM
impactor_particles.A2_vel -= A1_vel_COM

T = woma.A1_T_u_rho(target_particles.A1_u, target_particles.A1_rho, target_particles.A1_mat_id)
log_P = np.log10(target_particles.A1_P)
x, y, z = pos[:, 0], pos[:, 1], pos[:, 2]
r = np.sqrt(x ** 2 + y ** 2 + z ** 2)

target_particles.A1_mat_id[target_particles.A1_mat_id == 304] = 901

plt.scatter(T, log_P, s=5)
plt.xscale('log')
plt.ylim([12, 8])
plt.ylabel('Pressure ($\log$[Pa])')
plt.xlabel('Temperature (K)')
plt.savefig('PT sub neptune.png')
plt.close()

plt.scatter(r / R_earth, log_P, s=5)
plt.ylabel('Pressure ($\log$[Pa])')
plt.xlabel('r')
plt.savefig('Pr sub neptune.png')

#print(np.sum(target_particles.A1_m))
#print(np.sum(impactor_particles.A1_m))

print(np.max(sub_neptune.A1_rho))
print(np.max(sub_neptune.A1_P))

# import h5py
# with h5py.File(f'moon_sub_neptune_impact.hdf5', "w") as f:
#     woma.save_particle_data(
#         f,
#         np.append(target_particles.A2_pos, impactor_particles.A2_pos, axis=0),
#         np.append(target_particles.A2_vel, impactor_particles.A2_vel, axis=0),
#         np.append(target_particles.A1_m, impactor_particles.A1_m),
#         np.append(target_particles.A1_h, impactor_particles.A1_h),
#         np.append(target_particles.A1_rho, impactor_particles.A1_rho),
#         np.append(target_particles.A1_P, impactor_particles.A1_P),
#         np.append(target_particles.A1_u, impactor_particles.A1_u),
#         np.append(target_particles.A1_mat_id, impactor_particles.A1_mat_id),
#         boxsize=100 * R_earth,
#         file_to_SI=woma.Conversions(M_earth, R_earth, 1),
#     )

