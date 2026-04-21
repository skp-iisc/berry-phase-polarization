import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("alp/alat.txt", skiprows=1)
a, en = data[:, 0], data[:, 1]
plt.plot(a, en, 'o-')
plt.xlabel("Lattice parameter in Bohr")
plt.ylabel("Energy in Ry")
plt.tight_layout()
plt.savefig("alp/alat.png", dpi=300)
print("Check alp/alat.png.")

alat = 10.28        # CHECK from alp/alat.png minima
vol = alat**3
var_eps = 0.001     # electric field

# x components
F0, F1 = 0.00001735, 0.00315771
D0, D1 = -1.822773629679032E-005, 0.561041084984067

Z_star11 = (F1-F0)/(var_eps*np.sqrt(2))
print(f"Z* = {Z_star11}")

eps_inf = 1 + (4*np.pi/vol) * (D1-D0)/var_eps
print(f"eps_inf = {eps_inf}.")

Dr0, Dr1 = -1.822773629679032E-005, 0.502898501609083
Dtot0 = 188.995500475540 + Dr0
Dtot1 = 189.219245437475 + Dr1

eps_stat = 1 + (4*np.pi/vol) * (Dtot1-Dtot0)/var_eps
print(f"eps_stat = {eps_stat}.")
