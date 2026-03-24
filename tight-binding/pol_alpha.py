from ham_k import *

def calc_dyn_pol_alpha(T: float, t_max: float, N_k: int):
    '''Calculates dynamical polarization.

    Inputs
    ---

    T: time, t_max: maximum time, N_k: number of k-points.

    Returns
    ---

    time_array, P_dynamical.'''

    dt = 0.005          # Time step
    dk = 2 * np.pi / N_k
    
    D = np.array([np.exp(1j * dk / 3), 1.0, np.exp(-1j * dk / 3)])

    # 1
    H_init = fn_H_k_alpha(0.0, N_k)
    vals, vecs = np.linalg.eigh(H_init)
    v_k_dyn = vecs[:, :, 0]           # ground state wavefunction at t=0

    I_stack = np.tile(np.eye(3, dtype=complex), (N_k, 1, 1))
    
    # Time arrays
    times = np.arange(0, t_max + dt*0.1, dt)
    save_interval = 20 # input

    save_times = []
    P_dyn_wrapped = []

    # 2
    for step, t in enumerate(times):
        # Sliding parameter
        if t <= T:
            alpha = 2 * np.pi * (np.sin(np.pi * t / (2 * T)))**2
        else:
            alpha = 2 * np.pi
            
        if step % save_interval == 0:
            save_times.append(t)
            
            v_k_next = np.roll(v_k_dyn, shift=-1, axis=0)
            overlaps = np.sum(v_k_dyn.conj() * v_k_next * D, axis=1)
            P_dyn_wrapped.append(np.angle(np.prod(overlaps)) / (2 * np.pi))

        # 3
        H_all = fn_H_k_alpha(alpha, N_k)
        A = I_stack - 1j * (dt / 2) * H_all
        B = I_stack + 1j * (dt / 2) * H_all

        x = np.linalg.solve(B, v_k_dyn[..., np.newaxis])
        v_k_dyn = (A @ x)[..., 0]

    P_dyn = np.array(P_dyn_wrapped)
    P_dyn = np.unwrap(P_dyn * 2 * np.pi) / (2 * np.pi)

    P_dyn = P_dyn[0] - P_dyn
    
    return np.array(save_times), P_dyn

def calc_ad_pol_alpha(T: float, t_max: float, N_k: int):
    '''Calculates adiabatic polarization.

    Inputs
    ---

    T: time, t_max: maximum time, N_k: number of k-points.

    Returns
    ---

    time_array, P_adiabatic.'''

    dt = 0.005          # Time step
    dk = 2 * np.pi / N_k
    
    D = np.array([np.exp(1j * dk / 3), 1.0, np.exp(-1j * dk / 3)])
    
    # Time arrays
    times = np.arange(0, t_max + dt*0.1, dt)
    save_interval = 20

    save_times = []
    P_ad_wrapped = []

    # 2
    for step, t in enumerate(times):
        # Sliding parameter
        if t <= T:
            alpha = 2 * np.pi * (np.sin(np.pi * t / (2 * T)))**2
        else:
            alpha = 2 * np.pi
            
        if step % save_interval == 0:
            save_times.append(t)

            H_all = fn_H_k_alpha(alpha, N_k)
            vals, vecs = np.linalg.eigh(H_all)
            v_k_ad = vecs[:, :, 0]
            v_k_next_ad = np.roll(v_k_ad, shift=-1, axis=0)
            overlaps_ad = np.sum(v_k_ad.conj() * v_k_next_ad * D, axis=1)
            P_ad_wrapped.append(np.angle(np.prod(overlaps_ad)) / (2 * np.pi))

    P_ad = np.array(P_ad_wrapped)
    P_ad = np.unwrap(P_ad * 2 * np.pi) / (2 * np.pi)

    P_ad = P_ad[0] - P_ad
    
    return np.array(save_times), P_ad

