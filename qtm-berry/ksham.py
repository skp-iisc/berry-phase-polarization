# TODO: Refactor for spin-orbit calculatioons
from __future__ import annotations

from qtm.msg_format import shape_mismatch_msg

__all__ = ["KSHam"]
from collections.abc import Sequence
import numpy as np

from qtm.containers import FieldRType, get_WavefunG
from qtm.gspace import GkSpace
from qtm.pseudo import NonlocGenerator

from qtm.logger import qtmlogger
from qtm.msg_format import *


class KSHam:
    def __init__(
        self,
        gkspc: GkSpace,
        is_noncolin: bool,
        vloc: FieldRType,
        l_nloc: Sequence[NonlocGenerator],
        **kwargs
    ):
        if not isinstance(gkspc, GkSpace):
            raise TypeError(type_mismatch_msg("gkspc", gkspc, GkSpace))
        self.gkspc: GkSpace = gkspc
        self.ke_gk = get_WavefunG(gkspc, 1)((0.5 * self.gkspc.gk_norm2).astype("c16"))

        if not isinstance(is_noncolin, bool):
            raise TypeError(type_mismatch_msg("is_noncolin", is_noncolin, bool))
        self.is_noncolin: bool = is_noncolin

        if not isinstance(vloc, FieldRType):
            raise TypeError(type_mismatch_msg("vloc", vloc, FieldRType))
        if vloc.gspc is not gkspc.gwfn:
            raise ValueError(
                obj_mismatch_msg("vloc.gspc", vloc.gspc, "gkspc.gwfn", gkspc.gwfn)
            )

        if vloc.shape != (1 + is_noncolin,):
            if not is_noncolin and vloc.shape == ():
                pass
            else:
                raise ValueError(
                    value_mismatch_msg(
                        "vloc.shape",
                        vloc.shape,
                        f"{(1 + is_noncolin, )} when is_noncolin = {is_noncolin}",
                    )
                )
        self.vloc = vloc

        if not isinstance(l_nloc, Sequence) or any(
            not isinstance(nloc, NonlocGenerator) for nloc in l_nloc
        ):
            raise TypeError(type_mismatch_msg("l_nloc", l_nloc, NonlocGenerator))
        self.l_vkb_dij = []
        self.vnl_diag = 0
        for nloc in l_nloc:
            vkb, dij, vkb_diag = nloc.gen_vkb_dij(self.gkspc)
            self.l_vkb_dij.append((vkb, dij))
            self.vnl_diag += vkb_diag

        self.P_data = None
        self.Q_data = None
        
        efield_cart = kwargs.get('efield_cart', None)
        full_wfn = kwargs.get('full_wfn', None)
        
        if efield_cart is not None and full_wfn is not None:
            self._init_field_operator(
                efield_cart, full_wfn, kwargs['kpts'], 
                kwargs['kgrid_shape'], kwargs['ik'], 
                kwargs['crystal'], kwargs['current_kswfn']
            )

    @qtmlogger.time("KSHam:h_psi")
    def h_psi(self, l_psi: get_WavefunG, l_hpsi: get_WavefunG):
        # l_hpsi[:] = self.ke_gk * l_psi
        assert l_psi.shape == l_hpsi.shape, shape_mismatch_msg(
            "l_psi", "l_hpsi", l_psi, l_hpsi
        )
        l_psi = l_psi.reshape(-1)
        l_hpsi = l_hpsi.reshape(-1)

        np.multiply(self.ke_gk, l_psi, out=l_hpsi)
        for psi, hpsi in zip(l_psi, l_hpsi):
            psi_r = psi.to_r()
            psi_r *= self.vloc.data.ravel()
            hpsi += psi_r.to_g()

        for vkb, dij in self.l_vkb_dij:
            # proj = vkb.vdot(l_psi)
            # l_hpsi += (dij @ proj).T @ vkb

            proj = vkb.vdot(l_psi)
            proj = dij @ proj
            l_hpsi.zgemm(vkb.data.T, proj.T, 0, 1, 1.0, l_hpsi.data.T, 1.0)

        # Finite Field Modifier (PRL 89)
        if self.P_data is not None and self.Q_data is not None:
            proj_P = self.P_data.conj() @ l_psi.data.T
            proj_Q = self.Q_data.conj() @ l_psi.data.T
            l_hpsi.data[:] += (proj_P.T @ self.Q_data + proj_Q.T @ self.P_data)
            # l_hpsi.data[:] += proj_P.T @ self.Q_data

    def _init_field_operator(self, efield_cart, full_wfn, kpts, kgrid_shape, ik, crystal, current_kswfn):
        # Safety check: full_wfn should contain all k-points
        if full_wfn is None or len(full_wfn) != kpts.numkpts:
            self.Q_data = None
            return
            
        all_k = kpts.k_cryst
        k_cryst = all_k[:, ik]
        efield_frac = np.dot(crystal.reallat.latvec.T, efield_cart)
        
        numbnd = current_kswfn.numbnd
        nG_current = current_kswfn.evc_gk.data.shape[1]
        G_current = self.gkspc.g_cryst
        
        self.P_data = current_kswfn.evc_gk.data.copy()
        Q_data = np.zeros_like(self.P_data, dtype=np.complex128)
        
        def map_neighbor(wfn_neigh, u_list):
            mapped = np.zeros((numbnd, nG_current), dtype=np.complex128)
            G_neigh = wfn_neigh.gkspc.g_cryst
            neigh_dict = {tuple(G_neigh[:, i]): i for i in range(G_neigh.shape[1])}
            u_vec = np.array(u_list)
            for i_curr in range(nG_current):
                g_tar = tuple(G_current[:, i_curr] + u_vec)
                if g_tar in neigh_dict:
                    mapped[:, i_curr] = wfn_neigh.evc_gk.data[:, neigh_dict[g_tar]]
            return mapped

        for d in range(3):
            if abs(efield_frac[d]) < 1e-10: continue
            
            step = np.zeros(3)
            step[d] = 1.0 / kgrid_shape[d]      # Delta k
            
            k_target_p = (k_cryst + step) % 1.0
            k_target_m = (k_cryst - step) % 1.0
            
            ik_p = ik_m = None
            for i in range(kpts.numkpts):
                if np.allclose(all_k[:, i] % 1.0, k_target_p, atol=1e-5): ik_p = i
                if np.allclose(all_k[:, i] % 1.0, k_target_m, atol=1e-5): ik_m = i
                
            if ik_p is None or ik_m is None: continue
            
            wfn_p, wfn_m = full_wfn[ik_p], full_wfn[ik_m]
            umklapp_p = np.round(k_cryst - all_k[:, ik_p] + step).astype(int).tolist()
            umklapp_m = np.round(k_cryst - all_k[:, ik_m] - step).astype(int).tolist()
            
            bands = range(numbnd)
            S_p = current_kswfn.overlap(wfn_p, bands, bands, umklapp_p)
            S_m = current_kswfn.overlap(wfn_m, bands, bands, umklapp_m)
            
            S_p_inv = np.linalg.inv(S_p)
            S_m_inv = np.linalg.inv(S_m)
            
            mapped_wfn_p = map_neighbor(wfn_p, umklapp_p)
            mapped_wfn_m = map_neighbor(wfn_m, umklapp_m)
            
            Q_p = S_p_inv.T @ mapped_wfn_p
            Q_m = S_m_inv.T @ mapped_wfn_m
            
            prefactor = 1.0j * efield_frac[d] / (np.pi * step[d])
            Q_data += prefactor * (Q_p - Q_m)
            
        self.Q_data = Q_data
