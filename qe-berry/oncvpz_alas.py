import numpy as np
import matplotlib.pyplot as plt

print("Plotting lattice parameter vs. energy.")
data = np.loadtxt("alas/alat_pz.txt", skiprows=1)
a, en = data[:, 0], data[:, 1]
plt.plot(a, en, 'o-')
plt.xlabel("Lattice parameter in Bohr")
plt.ylabel("Energy in Ry")
plt.tight_layout()
# plt.savefig("alas/alas_alat_oncvpz.png", dpi=300)
plt.show()
# np.savetxt('alas/alat_oncvpz.txt', data, fmt='%.7f',)
# print("Check alas/alas_alat_oncvpz.png.")

alat = 10.58        # CHECK from alas/alas_alat.png minima
vol = alat**3
var_eps = 0.001

# x components
F0, F1 = -0.00000132, 0.00333573
D0, D1 = -6.495168322808666E-008, 0.749534388965161

Z_star11 = (F1-F0)/(var_eps*np.sqrt(2))
print(f"Z* = {Z_star11}")

eps_inf = 1 + (4*np.pi/vol) * (D1-D0)/var_eps
print(f"eps_inf = {eps_inf}.")

Dr0, Dr1 = -6.495168322808666E-008, 0.698764239341539
Dtot0 = 194.510933368796 + Dr0
Dtot1 = 194.748674058741 + Dr1

eps_stat = 1 + (4*np.pi/vol) * (Dtot1-Dtot0)/var_eps
print(f"eps_stat = {eps_stat}.")

