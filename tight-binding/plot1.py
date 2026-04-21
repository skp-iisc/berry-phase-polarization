from ham_k import *
import matplotlib.pyplot as plt

k_vals = np.linspace(-np.pi, np.pi, 200)
Ham_k = fn_H2_k_alpha(alpha=0, N_k=200)
eigvals = np.linalg.eigvalsh(Ham_k)
# energies[i, :] = np.sort(eigvals)
# print(Ham_k.shape, eigvals.shape)

plt.figure(figsize=(7, 5))
for ni in range(2):
    plt.plot(k_vals, eigvals[:, ni], 'k-')

plt.xlabel('$k$', fontsize=16)
plt.ylabel('Energy', fontsize=16)
plt.xlim(-np.pi, np.pi)
plt.ylim(-2.5, 2.5)
plt.xticks([-np.pi, 0, np.pi], [r'$-\pi$', '0', r'$\pi$'], fontsize=14)
plt.yticks([-2, -1, 0, 1, 2], fontsize=14)
plt.tick_params(direction='in', top=True, right=True, labelsize=16)
plt.tight_layout()
plt.savefig('figs/prb69fig1a.png', dpi=300)
plt.show()

Ham_k = fn_H_k_alpha(alpha=0, N_k=200)
eigvals = np.linalg.eigvalsh(Ham_k)

plt.figure(figsize=(7, 5))
for ni in range(3):
    plt.plot(k_vals, eigvals[:, ni], 'k-')

plt.xlabel('$k$', fontsize=16)
plt.ylabel('Energy', fontsize=16)
plt.xlim(-np.pi, np.pi)
plt.ylim(-2.5, 2.5)
plt.xticks([-np.pi, 0, np.pi], [r'$-\pi$', '0', r'$\pi$'], fontsize=14)
plt.yticks([-2, -1, 0, 1, 2], fontsize=14)
plt.tick_params(direction='in', top=True, right=True, labelsize=16)
plt.tight_layout()
plt.savefig('figs/prb69fig1.png', dpi=300)
plt.show()
