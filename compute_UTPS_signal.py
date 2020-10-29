from sys import argv
from copy import deepcopy
from itertools import product
import itertools

import pandas as pd
import numpy as np
import scipy as sp
import scipy.constants as spconst
from scipy.integrate import quad
from scipy.integrate import trapz
from matplotlib import pyplot as plt

DIPOLES_FILE = "globus\\GW_transition_dipole_files\\GW_transition_dipoles.csv"

# Functions for converting a dataframe to matrix in Hilbert or Liouville Space

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
        mu_mmp = mu_df.loc[(mu_label, *m_prime, *m)].to_numpy()[0]
    except KeyError as e:
        mu_mmp = 0
    
    return mu_nnp*delta(m, m_prime) - mu_mmp*delta(n, n_prime)

def make_mu_liouville_matrix(mu_df, rho_ordered):
    mus = {}
    for mu_label in ["x", "y", "z"]:
        mus[mu_label] = np.zeros(2*[len(rho_ordered)])
        for i, rho1 in enumerate(rho_ordered):
            for j, rho2 in enumerate(rho_ordered):
                mus[mu_label][i, j] = mu_liouville_element(mu_df, mu_label, rho1, rho2)
                                   
    return mus

def make_L_mol_liouville_matrix(energy_df, rho_ordered):
    L_mol = np.diag([energy_df.loc[rho[0]].to_numpy()[0] - energy_df.loc[rho[1]].to_numpy()[0] 
                                   for rho in rho_ordered])
    return L_mol

def mu_element(mu_df, mu_label, rho):
    n, m = rho
    try:
        mu_elm = mu_df.loc[(mu_label, *n, *m)].to_numpy()[0]
    except KeyError as e:
        mu_elm = 0
    
    return mu_elm

def make_mu_matrix(mu_df, state_order):
    mus = []
    for mu_label in ["x", "y", "z"]:
        mus.append(np.array([[mu_element(mu_df, mu_label, (elm1, elm2))
                                   for elm1 in state_order] for elm2 in state_order]))
    return np.array(mus)

def make_H_mol_matrix(energy_df, state_order):
    H_mol = np.zeros(2*[len(state_order)])
    for i, s1 in enumerate(state_order):
        H_mol[i,i] = energy_df.loc[s1].to_numpy()[0]
    return H_mol

# Functions to convert a Hilbert Space representation to Liouville and back

def hilbert_to_liouville_vector(hilbert_matrix):
    order_states = list(range(0, hilbert_matrix.shape[0]))
    order_pops = list(zip(order_states, order_states))
    order_cohrs = [(state1, state2) for state1, state2 in product(order_states, order_states)
                  if (state1, state2) not in order_pops]
    order = order_pops + order_cohrs
    
    liouville_vector = np.array([hilbert_matrix[i, j] for i,j in order])
    return liouville_vector

def liouville_vec_to_hilbert_matrix(liouville_vector):
    no_states = int(np.sqrt(len(liouville_vector)))
    order_states = list(range(0, no_states))
    order_pops = list(zip(order_states, order_states))
    order_cohrs = [(state1, state2) for state1, state2 in product(order_states, order_states)
                  if (state1, state2) not in order_pops]
    order = order_pops + order_cohrs
    
    hilbert_matrix = np.zeros(2*[no_states])
    for index, indices in enumerate(order):
        hilbert_matrix[indices[0], indices[1]] = liouville_vector[index]
    return hilbert_matrix

def hilbert_to_liouville_matrix(hilbert_matrix):
    no_states = hilbert_matrix.shape[0]
    order_states = list(range(0, no_states))
    order_pops = list(zip(order_states, order_states))
    order_cohrs = [(state1, state2) for state1, state2 in product(order_states, order_states)
                  if (state1, state2) not in order_pops]
    order = order_pops + order_cohrs
    
    liouville_matrix = np.zeros(2*[len(order)])
    for index, ipair in enumerate(order):
        irho_matrix = np.zeros(2*[no_states])
        irho_matrix[ipair[0], ipair[1]] = 1
        for jndex, jpair in enumerate(order):
            jrho_matrix = np.zeros(2*[no_states])
            jrho_matrix[jpair[0], jpair[1]] = 1
            
            com = hilbert_matrix @ irho_matrix - irho_matrix @ hilbert_matrix
            liouville_matrix[index, jndex] = np.trace(jrho_matrix.T @ com)
    return liouville_matrix

# Functions to compute an estimate of the Decay matrix and the zeroth order time evolution operator in liouville space

# guess for the Decay matrix in Liouville Space
def make_guess_from_energy(H_mol):
    no_states = H_mol.shape[0]
    order_states = list(range(0, no_states))
    order_pops = list(zip(order_states, order_states))
    order_cohrs = [(state1, state2) for state1, state2 in product(order_states, order_states)
                  if (state1, state2) not in order_pops]
    
    liouville_matrix = np.zeros(2*[len(order_pops + order_cohrs)])
    for index, ipair in enumerate(order_pops):
        for jndex, jpair in enumerate(order_pops):
            if index > jndex:
                liouville_matrix[jndex, index] = abs(H_mol[ipair] - H_mol[jpair])/100
    
    for jndex, jpair in enumerate(order_pops):
        liouville_matrix[jndex, jndex] = -sum(liouville_matrix[:jndex, jndex])
    
    offset = len(order_pops)
    for index, ipair in enumerate(order_cohrs):
        num = abs(H_mol[ipair[0], ipair[0]] - H_mol[ipair[1], ipair[1]])/100
        liouville_matrix[index + offset, index + offset] = num
    
    return liouville_matrix

def liou_time_op(L_0, t, t0):
    return sp.linalg.expm(-1j*L_0*(t - t0))

def liou_inverse(U_0):
    return np.linalg.inv(U_0)

# Functions that estimate input gaussian efields

def gaussian_pulse(z, t, z0, t0, amp, decay, freq, k):
    Z = z - z0
    T = t - t0
    return np.real(amp*np.exp(-decay*(Z/c-T)**2)*np.exp(-1j*freq*T-1j*k*Z))

def gaussian_Efield_z_prop(z, t, z0, t0, amp, decay, freq, k, theta):
    amp = gaussian_pulse(z, t, z0, t0, amp, decay, freq, k)
    x = amp*np.cos(theta)
    y = amp*np.sin(theta)
    z = 0
    return np.array([x, y, z])

# Functions to compute the lab frame response from molecular frame

def rot_ave_gamma_element(gamma, i, j, k, l):
    """
    Calculated the pure electronic rotationally averaged second hyperpolarizability in
    the lab frame using the molecular frame second order hyperpolarizability
    
    i,j,k,l are the lab frame cartiesian indices of rotationally ave gamma
    gamma is a 4th rank tensor with 4 indices in the molecular frame 
    
    Kwak 2015 Rigorous theory of molecular orientational nonlinear optics A7
    """

    def delta(i, j):
        """
        Python implementation of a delta function
        """
        if i == j:
            return 1
        else:
            return 0

    coordinates = [0, 1, 2]
    
    #Calculate pure electronic case
    del_fact0 = 4*delta(i, j)*delta(k, l) - delta(i, k)*delta(j, l) - delta(i, l)*delta(j, k)
    del_fact1 = 4*delta(i, k)*delta(j, l) - delta(i, j)*delta(k, l) - delta(i, l)*delta(j, k)
    del_fact2 = 4*delta(i, l)*delta(j, k) - delta(i, k)*delta(j, l) - delta(i, j)*delta(k, l)
    
    accumulator = 0
    called_indices = []
    for J in coordinates:
        for K in coordinates:
            accumulator += del_fact0*gamma[J,J,K,K] + del_fact1*gamma[J,K,J,K] + del_fact2*gamma[J,K,K,J]
            called_indices.append((J,J,K,K))
            called_indices.append((J,K,J,K))
            called_indices.append((J,K,K,J))
    accumulator /= 30
    return accumulator

def rot_ave_gamma(gamma):
    temp_gamma = np.zeros(gamma.shape, dtype=np.complex)
    
    coordinates = [0, 1, 2]
    for i,j,k,l in itertools.product([0,1,2], repeat=4):
        temp_gamma[i,j,k,l] = rot_ave_gamma_element(gamma, i, j, k, l)
    return temp_gamma
        
def rot_ave_UR_3_element(UR_3, j, k, l):
    """
    Calculates the rotationally averaged perturbed 3rd order time evolution response function.
    Which is used to calculate 3rd order time evolution operator and gamma
    """
    raise NotImplementedError

# Functions to compute the nth order Liouville Space time evolution operator

def tdot1(a, b):
    temp = np.tensordot(a, b, axes=[[-1], [1]])
    temp = np.swapaxes(temp, -3, -2)
    return temp

def tdot2(V, U):
    temp = np.tensordot(V, U, axes=[[-1], [0]])
    return temp

def tdot3(Uinv, V):
    temp = np.tensordot(V, Uinv, axes=[[-2], [1]])
    return temp

def U0_liou(L_0, t, t0):
    return sp.linalg.expm(-1j*L_0*(t - t0))

def U0_and_inv_liou(L_0, t, t0):
    temp = sp.linalg.expm(-1j*L_0*(t - t0))
    return temp, np.linalg.inv(temp)
    """
    except Exception as e:
        print("t, t0 = ", t, t0)
        U0_bound = max([abs(np.min(np.real(temp))), abs(np.max(np.real(temp)))])
        plt.imshow(np.real(temp), cmap='twilight', vmin=-U0_bound, vmax=U0_bound)
        plt.title("$Real \mathcal{U}_{0}$", fontsize=22)
        plt.colorbar()
        plt.show()
        
        U0_bound = max([abs(np.min(np.imag(temp))), abs(np.max(np.imag(temp)))])
        plt.imshow(np.imag(temp), cmap='twilight', vmin=-U0_bound, vmax=U0_bound)
        plt.title("$Imag \mathcal{U}_{0}$", fontsize=22)
        plt.colorbar()
        plt.show()
        raise e
        """

def perturbed_time_evo_resp_liou_gen(mu_liou_xyz, L_0, order=3):
    assert(order >= 1)
    def perturbed_liou_time_evo(*args):
        """
        parameters:
            t0, tau1, tau2, ..., tau_order, t_final
        """
        i = 0
        t0 = args[i]
        i += 1
        tau_next = args[i]
        i += 1
        
        U0, U0_inv = U0_and_inv_liou(L_0, tau_next, t0)
        temp = tdot3(U0_inv, tdot2(mu_liou_xyz, U0))
        
        for _ in range(order-1):
            tau_next = args[i]
            i += 1
            U0, U0_inv = U0_and_inv_liou(L_0, tau_next, t0)
            temp_next = tdot3(U0_inv, tdot2(mu_liou_xyz, U0))
            temp = tdot1(temp, temp_next)
        
        t_final = args[i]
        U0 = U0_liou(L_0, t_final, t0)
        temp = tdot3(U0, temp)
        return temp*1j**order
    return perturbed_liou_time_evo

def main():
    #atomic units
    c = spconst.physical_constants["inverse fine-structure constant"][0]
    amp, sigma, freq, k = 1*10**-6, 10**3, 1*10**-2, 0.0608/c
    decay = 1/sigma**2
    spec_amp, chirp = 1, 3
    t0, z0 = 0, 0

    #read and convert DIPOLES_FILE to df
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
    display(energy_df)
    state_order = energy_df.index.to_list()
    display(state_order)

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
    
    # Actual computation
    mu_xyz = make_mu_matrix(mu_df, state_order)

    mu_vec_xyz = []
    for index in range(3):
        mu_vec_xyz.append(hilbert_to_liouville_vector(mu_xyz[index]))
    mu_vec_xyz = np.array(mu_vec_xyz)

    mu_liou_xyz = []
    for index in range(3):
        mu_liou_xyz.append(hilbert_to_liouville_matrix(mu_xyz[index, :, :]))
    mu_liou_xyz = np.array(mu_liou_xyz)
    
    H_mol = make_H_mol_matrix(energy_df, state_order)

    L_mol = hilbert_to_liouville_matrix(H_mol)
    
    #estimate of the decay operator and zeroth order Liouville
    Decay = make_guess_from_energy(H_mol)

    L_0 = L_mol + 1j*Decay
    
    #estimate of rho at t0
    rho_t0 = np.zeros(2*[len(state_order)])
    rho_t0[0,0] = 1
    rho_t0_vec = hilbert_to_liouville_vector(rho_t0)

    #computing the time evolution operator and gamma
    UR_3 = perturbed_time_evo_resp_liou_gen(mu_liou_xyz, L_0, order=3)

    def gamma_R_3(ti, t0, t1, t2, tf):
        temp = np.tensordot(mu_vec_xyz, UR_3(ti, t0, t1, t2, tf), axes=[[1], [-2]])
        temp = np.tensordot(temp, rho_t0_vec, axes=[[-1], [0]])
        return temp

    def rot_ave_gamma_R_3(ti, t0, t1, t2, tf):
        return rot_ave_gamma(gamma_R_3(ti, t0, t1, t2, tf))
    
    #preparing the UTPS voltage function
    Eprobe = lambda t, tau, z: gaussian_Efield_z_prop(z, t, 0, tau, amp, decay, freq, k, 0)

    Edrive = lambda t, z: gaussian_Efield_z_prop(z, t, 0, 0, amp, decay, freq, k, np.pi/4)

    E_pulses = lambda t, tau, z: Eprobe(t, tau, z) + Edrive(t, z)

    def tdot_fields(U, v):
        temp = np.tensordot(U, v, axes=[[-1], [0]])
        return temp

    def integrand(ti, t0, t1, t2, tf, tau, z):
        return 1
        temp = tdot_fields(rot_ave_gamma_R_3(ti, t0, t1, t2, tf), E_pulses(t2 - ti, tau, z))
        temp = tdot_fields(temp, E_pulses(t1 - ti, tau, z))
        temp = tdot_fields(temp, E_pulses(t0 - ti, tau, z))
        temp = tdot_fields(temp, np.array([0,1,0]))
        return temp

    def YPolarization_lab(t_final, t_init, tau, z, no_time_steps):
        ts = np.linspace(t_init, t_final, no_time_steps)

        integrand_arr = np.zeros(3*[len(ts)], dtype=np.complex)
        integrand_arr.shape
        indicies = list(range(len(ts)))
        for i0, i1, i2 in itertools.product(indicies, repeat=3):
            t0, t1, t2 = ts[i0], ts[i1], ts[i2]
            if t_final >= t0 >= t1 >= t2:
                integrand_arr[i0, i1, i2] = integrand(t_init, t0, t1, t2, t_final, tau, z)
        return trapz(trapz(trapz(integrand_arr, ts), ts), ts)

    def Vy_lab(t_init, t_final, tau, z_init, z_final, k, no_t_steps, no_z_steps):
        ts = np.linspace(t_init, t_final, no_t_steps)
        dt = ts[0] - ts[1]
        zs = np.linspace(z_init, z_final, no_z_steps)
        Ts, Zs = np.meshgrid([ts, zs])
        
        yPol_arr = np.zeros([no_t_steps, no_z_steps], dtype=np.complex)
        for t, it in enumerate(ts):
            for z, iz in enumerate(zs):
                yPol_arr[it, iz] = YPolarization_lab(t_final, t_init, tau, z, no_t_steps)
        
        ddt_Py = np.gradient(yPol_arr, varargs=[dt], axis=0)
        Ey = 1j*trapz(ddt_Py*np.exp(1j*k*Zs), zs, axis=1)
        Vy = trapz(abs(Ey)**2, ts)
        return Vy
    
    # final signal calculation
    
    t_initial = 0
    t_final = 1*10**6 # ~30 ps in atomic units
    t_steps = 51
    
    z_init = 0
    z_final = 1*10**5 # ~5 microns
    z_steps = 51

    tau_init = -1*10**5
    tau_final = 1*10**5
    tau_steps = 31
    taus = np.linspace(t_init, t_final, tau_steps)

    with open("output_file.csv", "w") as f:
        f.write("tau (au.), OKE signal (arb.)\n")

    for tau in taus:
        UTPS_signal = Vy_lab(t_initial, t_final, tau, z_init, z_final, k, t_steps, z_steps)
        with open("output_file.csv", "a") as f:
            f.write(f"{tau},{UTPS_signal}\n")

