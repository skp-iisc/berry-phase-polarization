import numpy as np

def fn_H_k_alpha(alpha, N_k):
    """Constructs the 3 band Hamiltonian H(k).
    Inputs: alpha, N_k: number of k-points.
    Returns: H_k as a matrix."""
    # System parameters
    t_hop = 1.0
    Delta = -1.0

    beta = np.array([-2*np.pi/3, 0, 2*np.pi/3])

    k_vals = np.linspace(-np.pi, np.pi, N_k, endpoint=False)

    epsilon = Delta * np.cos(alpha - beta)
    H = np.zeros((N_k, 3, 3), dtype=complex)
    
    H[:, 0, 0] = epsilon[0]
    H[:, 1, 1] = epsilon[1]
    H[:, 2, 2] = epsilon[2]
    
    H[:, 0, 1] = H[:, 1, 0] = t_hop
    H[:, 1, 2] = H[:, 2, 1] = t_hop

    H[:, 0, 2] = t_hop * np.exp(-1j * k_vals)
    H[:, 2, 0] = t_hop * np.exp(1j * k_vals)
    return H

