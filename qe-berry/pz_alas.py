import numpy as np
import matplotlib.pyplot as plt

print("Plotting lattice parameter vs. energy.")
data = np.loadtxt("alas/alat.txt", skiprows=1)
a, en = data[:, 0], data[:, 1]
plt.plot(a, en, 'o-')
plt.xlabel("Lattice parameter in Bohr")
plt.ylabel("Energy in Ry")
plt.tight_layout()
# plt.savefig("alas/alas_alat_pz.png", dpi=300)
plt.show()
# np.savetxt('alas/alat_pz.txt', data, fmt='%.7f',)
# print("Check alas/alas_alat_pz.png.")

alat = 10.58        # CHECK from alas/alas_alat.png minima
vol = alat**3
var_eps = 0.001

# x components
F0, F1 = -0.00000035, 0.00306206
D0, D1 = -1.206023168480921E-005, 0.681767560375280

Z_star11 = (F1-F0)/(var_eps*np.sqrt(2))
print(f"Z* = {Z_star11}")

eps_inf = 1 + (4*np.pi/vol) * (D1-D0)/var_eps
print(f"eps_inf = {eps_inf}.")

Dr0, Dr1 = -1.206023168480921E-005, 0.912040868112379
Dtot0 = 134.661415409166 + Dr0
Dtot1 = 134.588671885614 + Dr1

eps_stat = 1 + (4*np.pi/vol) * (Dtot1-Dtot0)/var_eps
print(f"eps_stat = {eps_stat}.")

