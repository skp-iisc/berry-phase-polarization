import numpy as np
import os
from qtm.constants import RYDBERG
from qtm.lattice import RealLattice
from qtm.crystal import BasisAtoms, Crystal
from qtm.pseudo import UPFv2Data
from qtm.kpts import gen_monkhorst_pack_grid
from qtm.gspace import GSpace
from qtm.mpi import QTMComm
from qtm.dft import DFTCommMod, scf

from qtm.io_utils.dft_printers import print_scf_status
from qtm.logger import qtmlogger

from materials import MATERIALS
# qtmconfig.fft_backend = "mkl_fft"

from qtm.config import MPI4PY_INSTALLED
if MPI4PY_INSTALLED:
    from mpi4py.MPI import COMM_WORLD
else:
    COMM_WORLD = None

comm_world = QTMComm(COMM_WORLD)

# Only k-pt parallelization:
dftcomm = DFTCommMod(comm_world, comm_world.size, 1)
# Only band parallelization:
# dftcomm = DFTCommMod(comm_world, 1, 1)

mpgrid_shape = (4, 4, 4)
mpgrid_shift = (True, True, True)

ecut_wfn = 40 * RYDBERG
ecut_rho = 4 * ecut_wfn

conv_thr = 1e-8 * RYDBERG
diago_thr_init = 1e-2 * RYDBERG

from cryst import sc_crystal

crystal = sc_crystal('AlAs', alat=10.58)     # INPUT

kpts = gen_monkhorst_pack_grid(crystal, mpgrid_shape, mpgrid_shift)

grho = GSpace(crystal.recilat, ecut_rho)
gwfn = grho

numbnd = crystal.numel // 2

out = scf(
    dftcomm,
    crystal,
    kpts,
    grho,
    gwfn,
    numbnd,
    is_spin=False,
    is_noncolin=False,
    symm_rho=True,
    rho_start=None,
    occ_typ="fixed",
    conv_thr=conv_thr,
    diago_thr_init=diago_thr_init,
    iter_printer=print_scf_status,
)

scf_converged, rho, l_wfn_kgrp, en = out
