from ham_k import *

def calc_static_pol_el(E_field, N=200, max_iter=300, get_v_k=False):
    """
    Computes the stationary (adiabatic) polarization for a static field
    using the self-consistent iterative diagonalization of T_k.

    Inputs
    ---

    E_field: electric field strength, N: number of k-points.
    max_iter: maximum iteration for self-consistent calculation.
    get_v_k: get the wave functions.

    Returns
    ---

    P_static: static polarization.
    Note: the zero-field polarization is subtracted out to focus on the field-induced change.'''
    """
    if E_field == 0.0:
        if get_v_k:
            return None, 0.0
        else:
            return 0.0
        
    dk = 2*np.pi/N
    x_l = np.array([-1/3, 0, 1/3])
    
    # Intra-cell phases for accurate overlaps
    D_plus = np.exp(-1j * dk * x_l)
    D_minus = np.exp(1j * dk * x_l)
    
    H0 = fn_H_k_alpha(alpha=0.0, N_k=N)
    vals, vecs = np.linalg.eigh(H0)
    v_k = vecs[:, :, 0]
    
    # Self-consistent iteration
    for itr in range(max_iter):
        v_k_plus = np.roll(v_k, -1, axis=0)
        v_k_minus = np.roll(v_k, 1, axis=0)
        
        S_plus = np.sum(v_k.conj() * D_plus * v_k_plus, axis=1)
        S_minus = np.sum(v_k.conj() * D_minus * v_k_minus, axis=1)
        
        v_tilde_plus = (D_plus * v_k_plus) / S_plus[:, None]
        v_tilde_minus = (D_minus * v_k_minus) / S_minus[:, None]
        
        # Construct w_k operator
        w_k = 1j * E_field * N / (4 * np.pi) * (
            np.einsum('ni,nj->nij', v_tilde_plus, v_k.conj()) - 
            np.einsum('ni,nj->nij', v_tilde_minus, v_k.conj())
        )       # Einstein's summation for efficient outer product
        T_k = H0 + w_k + w_k.conj().transpose(0, 2, 1)
        
        vals, vecs = np.linalg.eigh(T_k)
        v_k_new = vecs[:, :, 0]
        
        # Include geometric phase in wave functions to prevent discontinuities
        phase = np.imag(np.log(np.sum(v_k.conj() * v_k_new, axis=1)))
        v_k_new = v_k_new * np.exp(-1j * phase)[:, None]
        
        if np.max(np.abs(v_k_new - v_k)) < 1e-7:    # Tolerance for convergence
            v_k = v_k_new
            break
        v_k = v_k_new
        if itr == max_iter - 1:
            print(f"Convergence NOT achieved for Static polarization after {max_iter} iterations!!!")
        
    v_k_plus = np.roll(v_k, -1, axis=0)
    overlaps = np.sum(v_k.conj() * D_plus * v_k_plus, axis=1)
    
    pol_static = np.imag(np.log(np.prod(overlaps))) / (2*np.pi)

    if get_v_k:
        return v_k, pol_static
    else:
        return pol_static


def calc_dyn_pol_el(E_max, T_ramp, N=200, dt=0.005, t_max=120):
    """
    Evolves the wavefunction under the TDSE as the field ramps up.
    Returns times and unwrapped dynamic polarization arrays.

    Inputs
    ---

    E_max: maximum electric field strength, T_ramp: ramp time, 
    N: number of k-points, dt: time step, t_max: maximum time.

    Returns
    ---

    time_array, P_dynamic. 
    Note that the zero-field polarization is subtracted out to 
    focus on the field-induced change.'''
    """
    dk = 2 * np.pi / N
    x_l = np.array([-1/3, 0, 1/3])
    D_plus = np.exp(-1j * dk * x_l)
    D_minus = np.exp(1j * dk * x_l)
    
    H0 = fn_H_k_alpha(alpha=0.0, N_k=N)
    vals, vecs = np.linalg.eigh(H0)
    v_k = vecs[:, :, 0]
    
    I_stack = np.tile(np.eye(3, dtype=complex), (N, 1, 1))
    times = np.arange(0, t_max + dt/2, dt)
    ang, saved_times = [], []
    save_interval = 20
    
    for step, t in enumerate(times):
        # Linear turn on of electric field
        E_field = E_max * min(t / T_ramp, 1.0)

        if step % save_interval == 0:
            v_k_plus = np.roll(v_k, -1, axis=0)
            overlaps = np.sum(v_k.conj() * D_plus * v_k_plus, axis=1)

            ang.append(np.imag(np.log(np.prod(overlaps))))
            saved_times.append(t)
            
        # dual states and overlaps
        v_k_plus = np.roll(v_k, -1, axis=0)
        v_k_minus = np.roll(v_k, 1, axis=0)
        S_plus = np.sum(v_k.conj() * D_plus * v_k_plus, axis=1)
        S_minus = np.sum(v_k.conj() * D_minus * v_k_minus, axis=1)
        
        v_tilde_plus = (D_plus * v_k_plus) / S_plus[:, None]
        v_tilde_minus = (D_minus * v_k_minus) / S_minus[:, None]
        
        w_k = 1j * E_field * N / (4 * np.pi) * (
            np.einsum('ni,nj->nij', v_tilde_plus, v_k.conj()) - 
            np.einsum('ni,nj->nij', v_tilde_minus, v_k.conj())
        )       # Einstein's summation for efficient outer product
        
        T_k = H0 + w_k + w_k.conj().transpose(0, 2, 1)
        
        # Unitary evolution
        A = I_stack - 1j * (dt / 2) * T_k
        B = I_stack + 1j * (dt / 2) * T_k
        x = np.linalg.solve(B, v_k[..., np.newaxis])
        v_k = (A @ x)[..., 0]
        
    P_dyn = np.unwrap(np.array(ang)) / (2*np.pi)
    P_dyn -= P_dyn[0]
    return np.array(saved_times), P_dyn

