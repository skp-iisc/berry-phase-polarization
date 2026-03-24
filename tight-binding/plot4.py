from pol_electric import *
import matplotlib.pyplot as plt

# print(calc_static_pol_el(0.01, N=200)[1])

E_max = 0.025
print("Simulating adiabatic (static) polarizations...")
E_arr = np.linspace(0, E_max, 50)
P_static_arr = np.array([calc_static_pol_el(E, N=200) for E in E_arr])
P_static_arr -= P_static_arr[0]
print("Simulating dynamic polarizations...")
times_40, P_dyn_40 = calc_dyn_pol_el(E_max, T_ramp=40, N=200, t_max=120)
times_80, P_dyn_80 = calc_dyn_pol_el(E_max, T_ramp=80, N=200, t_max=120)

# Interpolate static grid onto time arrays for the Adiatbatic lines
E_arr_40 = E_max * np.minimum(times_40 / 40.0, 1.0)
P_ad_40 = np.interp(E_arr_40, E_arr, P_static_arr)

E_arr_80 = E_max * np.minimum(times_80 / 80.0, 1.0)
P_ad_80 = np.interp(E_arr_80, E_arr, P_static_arr)

print("Simulating 100 k-point inset comparison...")
times_40_100k, P_dyn_40_100k = calc_dyn_pol_el(E_max, T_ramp=40, N=100, t_max=120)

print("Simulation done.")

plt.figure(figsize=(8, 6))
plt.plot(times_40, P_dyn_40, '-', lw=1, label='Dynamic (T=40)')
plt.plot(times_40, P_ad_40, 'k--', label='Adiabatic')
plt.plot(times_80, P_dyn_80, '-', lw=1, alpha=0.8, label='Dynamic (T=80)')
plt.plot(times_80, P_ad_80, 'k--', alpha=0.5)

plt.xlabel('Time', fontsize=12)
plt.ylabel('Polarization', fontsize=12)
plt.xlim(0, 120)
plt.ylim(0, 0.0025)
plt.xticks([0, 40, 80, 120], fontsize=10)
plt.yticks([0, 0.001, 0.002], fontsize=10)
plt.legend(loc='upper left', fontsize=10)
plt.tick_params(direction='in', top=True, right=True)
plt.tight_layout()
plt.savefig('figs/prb69fig4a.png', dpi=300)
plt.show()
print("Polarization evolution in ramped electric field is plotted.")

plt.figure(figsize=(8, 6))
mask = (times_40 >= 38) & (times_40 <= 75)
mask_100 = (times_40_100k >= 38) & (times_40_100k <= 75)

plt.plot(times_40[mask], P_dyn_40[mask], '-', label='Dynamic (200 k points)')
plt.plot(times_40_100k[mask_100], P_dyn_40_100k[mask_100], color='gray', linestyle='-.', label='Dynamic (100 k points)')
plt.plot(times_40[mask], P_ad_40[mask], 'k--', label='Adiabatic')

plt.xlabel('Time', fontsize=12)
plt.ylabel('Polarization', fontsize=12)
plt.xlim(38, 75)
plt.ylim(0.00215, 0.00230)
plt.xticks([40, 60], fontsize=10)
plt.yticks([0.00215, 0.00220, 0.00225], fontsize=10)
plt.legend(loc='lower right', fontsize=10)
plt.tick_params(direction='in', top=True, right=True)
plt.tight_layout()
plt.savefig('figs/prb69fig4b.png', dpi=300)
plt.show()
print("Remnant oscillations (zoomed-in) is plotted.")

