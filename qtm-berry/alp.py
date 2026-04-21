import numpy as np
import os
import itertools
from qtm.constants import RYDBERG
from qtm.lattice import RealLattice
from qtm.crystal import BasisAtoms, Crystal
from qtm.pseudo import UPFv2Data
from qtm.kpts import gen_monkhorst_pack_grid, KList
from qtm.gspace import GSpace
from qtm.mpi import QTMComm
from qtm.dft import DFTCommMod, scf

from qtm.io_utils.dft_printers import print_scf_status
from qtm.logger import qtmlogger

from materials import MATERIALS
from cryst import sc_crystal

from qtm.config import MPI4PY_INSTALLED
if MPI4PY_INSTALLED:
    from mpi4py.MPI import COMM_WORLD
else:
    COMM_WORLD = None

comm_world = QTMComm(COMM_WORLD)
dftcomm = DFTCommMod(comm_world, comm_world.size, 1)

mpgrid_shape = (4, 4, 4)
ecut_wfn = 40 * RYDBERG
ecut_rho = 4 * ecut_wfn

conv_thr = 1e-7 * RYDBERG
diago_thr_init = 1e-2 * RYDBERG

crystal = sc_crystal('AlP', alat=10.28)

kpts = gen_monkhorst_pack_grid(crystal, mpgrid_shape, (True, True, True),
                                   use_symm=False, is_time_reversal=False)

grho = GSpace(crystal.recilat, ecut_rho)
gwfn = grho
numbnd = (crystal.numel//2) #+ 4     # 4 conduction bands

# Define field magnitude strictly in Hartree atomic units
E0_ha = 0.001

if comm_world.rank == 0:
    print(f"\n[PHASE 1: Extracting Exact P(0) Baseline]")

out_zero = scf(
    dftcomm, crystal, kpts, grho, gwfn, numbnd,
    is_spin=False, is_noncolin=False, 
    symm_rho=True,  # Protect the cubic ground state symmetry
    rho_start=None, occ_typ="fixed", 
    conv_thr=conv_thr, diago_thr_init=diago_thr_init,
    efield_cart=[0.0, 0.0, 0.0], 
    kgrid_shape=mpgrid_shape,    
    iter_printer=print_scf_status,
)

rho_0 = out_zero[1]
P_zero_vec = out_zero[4]  # The added P_el_cart return value

if comm_world.rank == 0:
    print(f"\n[PHASE 2: Applying X-Directed E-Field = {E0_ha} Ha a.u.]")

# conv_thr = 1e-5 * RYDBERG

out_field = scf(
    dftcomm, crystal, kpts, grho, gwfn, numbnd,
    is_spin=False, is_noncolin=False, 
    symm_rho=False,  # Field breaks cubic symmetry   
    rho_start=rho_0, # Restart from unperturbed density for faster convergence   
    occ_typ="fixed", 
    conv_thr=conv_thr, diago_thr_init=diago_thr_init,
    efield_cart=[E0_ha, 0.0, 0.0], 
    kgrid_shape=mpgrid_shape,
    iter_printer=print_scf_status,
    maxiter=10
)

scf_converged_E = out_field[0]
P_E_cart = out_field[4]

if comm_world.rank == 0:
    print(f"\nSCF Converged under field: {scf_converged_E}")
    print("========================================")
    print(f"P(0,0,0) Baseline      : {P_zero_vec} e/bohr^2")
    print(f"P(E,0,0) Perturbed     : {P_E_cart} e/bohr^2")
    print("----------------------------------------")
    
    # along the x-axis
    delta_P_x = P_E_cart[0] - P_zero_vec[0]
    chi_xx = delta_P_x / E0_ha

    eps_inf = 1.0 + 4.0 * np.pi * chi_xx
    
    print(f"Delta P_x (Linear): {delta_P_x:.8f}.")
    print(f"Dielectric eps_inf: {eps_inf:.4f}.")
    print("========================================")
