import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("alas/alat.txt", skiprows=1)
a, en = data[:, 0], data[:, 1]
plt.plot(a, en, 'o-')
plt.xlabel("Lattice parameter in Bohr")
plt.ylabel("Energy in Ry")
plt.tight_layout()
plt.savefig("alas/alat.png", dpi=300)
print("Check alas/alat.png.")

alat = 10.58        # CHECK from alas/alat.png minima
vol = alat**3
var_eps = 0.001     # electric field

# x components
F0, F1 = -0.00003044, 0.00295427
D0, D1 = -1.631019853032853E-005, 0.699649852543584

Z_star11 = (F1-F0)/(var_eps*np.sqrt(2))
print(f"Z* = {Z_star11}")

eps_inf = 1 + (4*np.pi/vol) * (D1-D0)/var_eps
print(f"eps_inf = {eps_inf}.")

Dr0, Dr1 = -1.631019853032853E-005, 0.636438065692201
Dtot0 = 194.510933368796 + Dr0
Dtot1 = 194.723340248345 + Dr1

eps_stat = 1 + (4*np.pi/vol) * (Dtot1-Dtot0)/var_eps
print(f"eps_stat = {eps_stat}.")
