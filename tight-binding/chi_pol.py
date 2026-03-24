from pol_electric import *

def calc_chi_real_imag(E0, omegas, dE=0.0005, N=100, dt=0.01, t_max=250, delta=0.04):
    """
    Computes the real and imaginary parts of the dielectric function chi.

    Inputs
    ---

    E0: electric field strength, omega: frequency, dE, N, dt, t_max, delta.

    Returns
    ---

    re_chi, im_chi: real and imaginary parts of chi.
    """
    dk = 2 * np.pi / N
    x_l = np.array([-1/3, 0, 1/3])
    D_plus = np.exp(-1j * dk * x_l)
    D_minus = np.exp(1j * dk * x_l)
    H0 = fn_H_k_alpha(alpha=0.0, N_k=N)

    # 1
    v_k_initial, P_init = calc_static_pol_el(E0 + dE, N=N, get_v_k=True)
    P_static_E0 = calc_static_pol_el(E0, N=N)
    dP_dE = (P_init - P_static_E0) / dE

    I_stack = np.tile(np.eye(3, dtype=complex), (N, 1, 1))
    times = np.arange(0, t_max + dt/2, dt)
    delta_P_t = np.zeros(len(times))
    v_k = v_k_initial
    
    # 2
    for step, t in enumerate(times):
        v_k_plus = np.roll(v_k, -1, axis=0)
        overlaps = np.sum(v_k.conj() * D_plus * v_k_plus, axis=1)
        P_current = np.angle(np.prod(overlaps)) / (2 * np.pi)
        
        # Delta P(t) tracking
        diff = P_current - P_static_E0
        delta_P_t[step] = (diff + 0.5) % 1.0 - 0.5
        
        v_k_minus = np.roll(v_k, 1, axis=0)
        S_plus = np.sum(v_k.conj() * D_plus * v_k_plus, axis=1)
        S_minus = np.sum(v_k.conj() * D_minus * v_k_minus, axis=1)

        v_tilde_plus = (D_plus * v_k_plus) / S_plus[:, None]
        v_tilde_minus = (D_minus * v_k_minus) / S_minus[:, None]

        w_k = 1j * E0 * N / (4 * np.pi) * (
            np.einsum('ni,nj->nij', v_tilde_plus, v_k.conj()) -
            np.einsum('ni,nj->nij', v_tilde_minus, v_k.conj())
        )
        T_k = H0 + w_k + w_k.conj().transpose(0, 2, 1)

        A = I_stack - 1j * (dt / 2) * T_k
        B = I_stack + 1j * (dt / 2) * T_k
        x = np.linalg.solve(B, v_k[..., np.newaxis])
        v_k = (A @ x)[..., 0]

    # 3
    times_grid, omegas_grid = np.meshgrid(times, omegas)
    integrand = delta_P_t * np.exp((1j * omegas_grid - delta) * times_grid)
    dP_w = dt * (0.5*integrand[:, 0] + np.sum(integrand[:, 1:-1], axis=1) 
                 + 0.5*integrand[:, -1])
    
    # 4
    re_chi = dP_dE - (omegas / dE) * np.imag(dP_w)
    im_chi = (omegas / dE) * np.real(dP_w)
    
    # 5
    re_chi = 1.0 + 4 * np.pi * re_chi
    im_chi = 4 * np.pi * im_chi
    
    return re_chi, im_chi

