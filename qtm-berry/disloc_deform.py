from cryst import *

from qtm.constants import RYDBERG, FPI, TPI
from qtm.kpts import KList, gen_monkhorst_pack_grid
from qtm.gspace import GSpace
from qtm.mpi import QTMComm
from qtm.dft import DFTCommMod, DFTConfig, scf
from qtm.io_utils.dft_printers import print_scf_status
from qtm.config import MPI4PY_INSTALLED

if MPI4PY_INSTALLED:
    from mpi4py.MPI import COMM_WORLD
else:
    COMM_WORLD = None

comm_world = QTMComm(COMM_WORLD)

# Only k-pt parallelization:
dftcomm = DFTCommMod(comm_world, comm_world.size, 1)

from materials import MATERIALS

def setup_crystal(material_name, alat=None, cation_disp=None, anion_disp=None,
                  strain_tensor=None):
    """
    Set up crystal for III-V semiconductors.

    Inputs
    ----------
    material_name : str - 'AlAs', 'GaAs', 'GaP', or 'AlP'
    alat : float or None - lattice constant in Bohr
    cation_disp : (3,) array - cation displacement in crystal coords
    anion_disp : (3,) array - anion displacement in crystal coords
    strain_tensor : (3,3) array - strain eta (a_i --> (I + eta) a_i)

    Returns
    ----------
    Crystal(reallat, [cation_atoms, anion_atoms])
    """
    mat = MATERIALS[material_name]
    if alat is None:
        alat = mat['alat']

    # FCC primitive vectors (in units of alat)
    a1 = np.array([-0.5, 0.0, 0.5])
    a2 = np.array([0.0, 0.5, 0.5])
    a3 = np.array([-0.5, 0.5, 0.0])

    if strain_tensor is not None:
        I_plus_eta = np.eye(3) + np.array(strain_tensor)
        a1, a2, a3 = I_plus_eta @ a1, I_plus_eta @ a2, I_plus_eta @ a3

    reallat = RealLattice.from_alat(alat=alat, a1=a1, a2=a2, a3=a3)

    r_cation = np.array([0.0, 0.0, 0.0])
    r_anion = np.array([0.25, 0.25, 0.25])
    if cation_disp is not None:
        r_cation += np.array(cation_disp)
    if anion_disp is not None:
        r_anion += np.array(anion_disp)

    pp_dir = os.path.dirname(os.path.abspath(__file__))
    cation_pp = UPFv2Data.from_file(os.path.join(pp_dir, mat['cation_pp']))
    anion_pp = UPFv2Data.from_file(os.path.join(pp_dir, mat['anion_pp']))

    cation_atoms = BasisAtoms(mat['cation'], cation_pp, mat['cation_mass'],
                              reallat, r_cation.reshape(3, 1))
    anion_atoms = BasisAtoms(mat['anion'], anion_pp, mat['anion_mass'],
                             reallat, r_anion.reshape(3, 1))

    return Crystal(reallat, [cation_atoms, anion_atoms])

def modified_sc_crystal(material_name, alat=None, cation_disp=None, anion_disp=None,
                  strain_tensor=None):
    """
    Set up modified crystal for III-V semiconductors.

    Inputs
    ----------
    material_name : str - 'Si', 'AlAs', 'GaAs', 'GaP', or 'AlP'
    alat : float or None - lattice constant in Bohr

    Returns
    ----------
    Crystal(reallat, [cation_atoms, anion_atoms])
    """
    mat = MATERIALS[material_name]
    if alat is None:
        alat = mat['alat']

    # Simple Cubic (SC) primitive vectors (in units of alat)
    a1 = np.array([1.0, 0.0, 0.0])
    a2 = np.array([0.0, 1.0, 0.0])
    a3 = np.array([0.0, 0.0, 1.0])

    if strain_tensor is not None:
        I_plus_eta = np.eye(3) + np.array(strain_tensor)
        a1, a2, a3 = I_plus_eta @ a1, I_plus_eta @ a2, I_plus_eta @ a3

    reallat = RealLattice.from_alat(alat=alat, a1=a1, a2=a2, a3=a3)

    rc1 = np.array([-0.125, -0.125, -0.125])
    rc2 = np.array([0.375,  0.375, -0.125])
    rc3 = np.array([0.375, -0.125,  0.375])
    rc4 = np.array([-0.125,  0.375,  0.375])
    ra1 = np.array([0.125,  0.125,  0.125])
    ra2 = np.array([0.625,  0.625,  0.125])
    ra3 = np.array([0.625,  0.125,  0.625])
    ra4 = np.array([0.125,  0.625,  0.625])

    if cation_disp is not None:
        rc1 += np.array(cation_disp)
    if anion_disp is not None:
        ra1 += np.array(anion_disp)

    rc = np.column_stack([rc1, rc2, rc3, rc4])
    ra = np.column_stack([ra1, ra2, ra3, ra4])

    pp_dir = os.path.dirname(os.path.abspath(__file__))
    cation_pp = UPFv2Data.from_file(os.path.join(pp_dir, mat['cation_pp']))
    anion_pp = UPFv2Data.from_file(os.path.join(pp_dir, mat['anion_pp']))

    cations = BasisAtoms(mat['cation'], cation_pp, mat['cation_mass'],
                        reallat, rc)
    anions = BasisAtoms(mat['anion'], anion_pp, mat['anion_mass'],
                        reallat, ra)

    return Crystal(reallat, [cations, anions])

def run_scf(crystal, kgrid_shape=(4, 4, 4), ecut_ry=25, numbnd=None,
            conv_thr=1e-7):
    """
    Run SCF on FULL (unsymmetrized) k-grid needed for Berry phase.

    Returns: l_wfn_kgrp, kpts, rho, en.
    """
    comm_world = QTMComm(COMM_WORLD)
    dftcomm = DFTCommMod(comm_world, comm_world.size, 1)

    # Full k-grid: no symmetry, no time-reversal
    kpts = gen_monkhorst_pack_grid(
        crystal, kgrid_shape,
        shifts=(False, False, False), use_symm=False, is_time_reversal=False
    )

    ecut_wfn = ecut_ry * RYDBERG
    grho = GSpace(crystal.recilat, 4 * ecut_wfn)
    gwfn = grho
    if numbnd is None:
        numbnd = crystal.numel // 2

    dftconfig1 = DFTConfig()
    dftconfig1.eigsolve_method == "davidson"
    out = scf(
        dftcomm, crystal, kpts, grho, gwfn, numbnd,
        is_spin=False, is_noncolin=False, symm_rho=True,
        occ_typ='fixed', conv_thr=conv_thr * RYDBERG,
        diago_thr_init=1e-2 * RYDBERG,
        dftconfig=dftconfig1,
        iter_printer=print_scf_status,
        maxiter=100,
    )
    scf_converged, rho, l_wfn_kgrp, en = out
    if not scf_converged:
        print("WARNING: SCF did not converge!")
    return l_wfn_kgrp, kpts, rho, en

def sort_kpts_into_strings(kpts, grid_shape, direction):
    """
    Organize k-points into Berry phase strings along a reciprocal direction.

    Returns list of lists of k-point indices.
    """
    N = grid_shape
    nk = kpts.numkpts
    assert nk == N[0] * N[1] * N[2], \
        f"Expected {N[0]*N[1]*N[2]} k-points, got {nk}."

    all_k = kpts.k_cryst  # (3, nk)
    perp_dirs = [d for d in range(3) if d != direction]

    # Group by perpendicular components
    strings = {}
    for ik in range(nk):
        k = all_k[:, ik]
        perp_key = tuple(round(k[d] * N[d]) % N[d] for d in perp_dirs)
        if perp_key not in strings:
            strings[perp_key] = []
        strings[perp_key].append(ik)

    # Sort each string by parallel component
    result = []
    for key in sorted(strings.keys()):
        indices = strings[key]
        indices.sort(key=lambda ik: all_k[direction, ik] % 1.0)
        assert len(indices) == N[direction], \
            f"Expected {N[direction]} per string, got {len(indices)}"
        result.append(indices)

    n_perp = N[perp_dirs[0]] * N[perp_dirs[1]]
    assert len(result) == n_perp
    return result

def calc_berry_phase(l_wfn, kpts, crystal, grid_shape, direction):
    """
    Compute the string-averaged Berry phase along reciprocal direction i.

    phi^(i) = (1/N_perp) * sum_strings Im ln prod det S(k_j, k_{j+1})

    Returns phi in radians.
    """
    numbnd_occ = crystal.numel // 2
    strings = sort_kpts_into_strings(kpts, grid_shape, direction)
    
    Gam = 0.0
    for string in strings:
        product = 1.0 + 0.0j
        N_str = len(string)
        for j in range(N_str):
            ik1 = string[j]
            ik2 = string[(j+1) % N_str]

            # Umklapp for wrapping around BZ
            umklapp = [0, 0, 0]
            if j == N_str-1:
                umklapp[direction] = 1

            wfn1 = l_wfn[ik1][0]
            wfn2 = l_wfn[ik2][0]
            bands = range(numbnd_occ)
            S = wfn1.overlap(wfn2, bands, bands, umklapp)
            product *= np.linalg.det(S)
        
        Gam += np.imag(np.log(product))
    
    return Gam/len(strings)
