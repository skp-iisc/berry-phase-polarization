from disloc_deform import *

material = "AlAs"
grid_shape = (4, 4, 4)
ecut_ry = 40.0
delta = 0.005

"""
Born effective charge Z*.

Z* = Z_ion - (1/pi) * dGam / (2*delta)

where dGam is the Berry phase change when the cation is displaced
by +/- delta in crystal coordinates along direction 0.
"""

print(f"\n{'='*60}")
print(f"Computing Z* for {material}")
print(f"{'='*60}")

cryst_p = modified_sc_crystal(material, cation_disp=[delta, 0, 0])
wfn_p, kpts_p, _, _ = run_scf(cryst_p, grid_shape, ecut_ry)
Gam_p = np.array([calc_berry_phase(wfn_p, kpts_p, cryst_p, grid_shape, d) 
                    for d in range(3)])
cryst_m = modified_sc_crystal(material, cation_disp=[-delta, 0, 0])
wfn_m, kpts_m, _, _ = run_scf(cryst_m, grid_shape, ecut_ry)
Gam_m = np.array([calc_berry_phase(wfn_m, kpts_m, cryst_m, grid_shape, d) 
                    for d in range(3)])
dGam = Gam_p - Gam_m

Z_ion = cryst_p.l_atoms[0].ppdata.valence
Z_star = Z_ion - (1.0/np.pi) * dGam[0]/(2*delta)

print(f"delta = {delta} => del_Gamma = {dGam}.")
print(f"Z_ion = {Z_ion}, Z* = {Z_star}.")
