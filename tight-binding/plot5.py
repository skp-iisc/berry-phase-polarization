from pol_electric import *
import matplotlib.pyplot as plt

print("Simulating dynamic polarizations...")
E_max = 0.05
times_40_800, P_dyn_40_800 = calc_dyn_pol_el(E_max, T_ramp=40, N=800, t_max=120)
times_80_800, P_dyn_80_800 = calc_dyn_pol_el(E_max, T_ramp=80, N=800, t_max=120)
print("Simulating 400 k-point inset comparison...")
times_40_400, P_dyn_40_400 = calc_dyn_pol_el(E_max, T_ramp=40, N=400, t_max=120)
print("Simulation done.")

plt.figure(figsize=(8, 6))
plt.plot(times_40_800, P_dyn_40_800, '-', lw=1, label='Dynamic (T=40)')
plt.plot(times_80_800, P_dyn_80_800, '-', lw=1, alpha=0.8, label='Dynamic (T=80)')

plt.xlabel('Time', fontsize=12)
plt.ylabel('Polarization', fontsize=12)
plt.xlim(0, 120)
plt.ylim(0, 0.005)
plt.xticks([0, 40, 80, 120], fontsize=10)
plt.yticks([0, 0.001, 0.002, 0.003, 0.004], fontsize=10)
plt.legend(loc='upper left', fontsize=10)
plt.tick_params(direction='in', top=True, right=True)
plt.tight_layout()
plt.savefig('figs/prb69fig5a.png', dpi=300)
plt.show()
print("Polarization evolution in ramped electric field is plotted.")

plt.figure(figsize=(8, 6))
mask = (times_40_800 >= 38) & (times_40_800 <= 75)
mask_400 = (times_40_400 >= 38) & (times_40_400 <= 75)

plt.plot(times_40_800[mask], P_dyn_40_800[mask], '-', label='Dynamic (800 k points)')
plt.plot(times_40_400[mask_400], P_dyn_40_400[mask_400], '.', color='gray', ms=4, label='Dynamic (400 k points)')

plt.xlabel('Time', fontsize=12)
plt.ylabel('Polarization', fontsize=12)
plt.xlim(38, 75)
plt.ylim(0.0043, 0.0046)
plt.xticks([40, 60], fontsize=10)
plt.yticks([0.0043, 0.0044, 0.0045], fontsize=10)
plt.legend(loc='lower right', fontsize=10)
plt.tick_params(direction='in', top=True, right=True)
plt.tight_layout()
plt.savefig('figs/prb69fig5b.png', dpi=300)
plt.show()
print("Remnant oscillations (zoomed-in) is plotted.")

