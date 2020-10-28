from sys import argv
from copy import deepcopy
from itertools import product

import pandas as pd
import numpy as np
import scipy as sp

DIPOLES_FILE = "globus\\GW_transition_dipole_files\\GW_transition_dipoles.csv"

def delta(state1, state2):
    if all([elm1 == elm2 for elm1, elm2 in zip(state1, state2)]):
        return 1
    else:
        return 0

def mu_liouville_element(mu_df, mu_label, rho_1, rho_2):
    n_prime, m_prime = rho_1
    n, m = rho_2
    try:
        mu_nnp = mu_df.loc[(mu_label, *n, *n_prime)].to_numpy()[0]
    except KeyError as e:
        mu_nnp = 0
    try:
        mu_mmp = mu_df.loc[(mu_label, *m, *m_prime)].to_numpy()[0]
    except KeyError as e:
        mu_mmp = 0
    
    return mu_nnp*delta(m, m_prime) - mu_mmp*delta(n, n_prime)

def make_mu_liouville_matrix(mu_df, rho_ordered):
    mus = {}
    for mu_label in ["x", "y", "z"]:
        mus[mu_label] = np.array([[mu_liouville_element(mu_df, mu_label, rho1, rho2)
                                   for rho1 in rho_ordered] for rho2 in rho_ordered])
    return mus

def make_L_mol_liouville_matrix(energy_df, rho_ordered):
    L_mol = np.diag([energy_df.loc[rho[0]].to_numpy()[0] - energy_df.loc[rho[1]].to_numpy()[0] 
                                   for rho in rho_ordered])
    return L_mol

def main(input_dipoles_path):
    df = pd.read_csv(DIPOLES_FILE, index_col=0)
    mu_df2 = pd.read_csv(DIPOLES_FILE, index_col=[1,4,5,6,7,8,9])
    df = df.fillna(0)
    df = df.drop(columns=["mu_sym","mu_spin"])
    df['B_spin'] = df['B_spin'].astype(int)
    df['C_spin'] = df['C_spin'].astype(int)

    benergy_df = df.drop(columns=["mu_label", "C_energy", "C_state_no", "C_sym", "C_spin", "moment"])
    benergy_df = benergy_df.drop_duplicates()
    benergy_df = benergy_df.rename(columns={"B_state_no": "state_no", "B_sym": "sym", "B_spin": "spin", "B_energy": "energy"})

    cenergy_df = df.drop(columns=["mu_label", "B_energy", "B_state_no", "B_sym", "B_spin", "moment"])
    cenergy_df = cenergy_df.drop_duplicates()
    cenergy_df = cenergy_df.rename(columns={"C_state_no": "state_no", "C_sym": "sym", "C_spin": "spin", "C_energy": "energy"})

    energy_df = pd.concat([benergy_df, cenergy_df])
    energy_df = energy_df.drop_duplicates(subset=["state_no", "sym", "spin"])
    energy_df = energy_df.sort_values("energy")
    energy_df = energy_df.set_index(["state_no", "sym", "spin"])
    state_order = energy_df.index.to_list()

    rho_populations = list(zip(state_order, state_order))
    rho_coherences = [ (state1, state2) for state1, state2 in product(state_order, state_order)
                      if (state1, state2) not in rho_populations]
    rho_ordered = rho_populations + rho_coherences

    mu_df = df.drop(columns=["B_energy", "C_energy"])
    mu_df_adj = mu_df.rename(columns={"C_state_no": "B_state_no", "C_sym": "B_sym",
                                      "C_spin": "B_spin", "B_state_no": "C_state_no",
                                      "B_sym": "C_sym", "B_spin": "C_spin"})

    mu_df = pd.concat([mu_df, mu_df_adj])
    mu_df = mu_df.drop_duplicates(subset=["mu_label", "B_state_no", "B_sym", "B_spin", "C_state_no", "C_sym", "C_spin"])
    mu_df = mu_df.set_index(["mu_label", "B_state_no", "B_sym", "B_spin", "C_state_no", "C_sym", "C_spin"]).sort_index()
    
    mu_liou_xyz = make_mu_liouville_matrix(mu_df, rho_ordered)
    L_mol = make_L_mol_liouville_matrix(energy_df, rho_ordered)
    
if __name__ == "__main__":
    main()
