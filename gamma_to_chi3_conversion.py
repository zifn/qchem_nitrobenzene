import os
import itertools
import pickle
import sympy as sp
import numpy as np
import pandas as pd
from sympy import init_printing
import scipy.constants as constants
init_printing()

def euler_integration(R):
    return sp.integrate(sp.integrate(sp.integrate(R, (alph, -sp.pi, sp.pi)), (beta, -sp.pi/2, sp.pi/2)), (gamma, -sp.pi, sp.pi))

def delta(i, j):
    """
    Python implementation of a delta function
    """
    if i == j:
        return 1
    else:
        return 0
    
def rot_ave_gamma(gamma, i, j, k, l):
    """
    Calculated the pure electronic rotationally averaged second hyperpolarizability in
    the lab frame using the molecular frame second order hyperpolarizability
    
    i,j,k,l are the lab frame cartiesian indices of rotationally ave gamma
    gamma is a 4th rank tensor with 4 indices in the molecular frame 
    
    Kwak 2015 Rigorous theory of molecular orientational nonlinear optics A7
    """
    coordinates = [0, 1, 2]
    
    #Calculate pure electronic case
    del_fact0 = 4*delta(i, j)*delta(k, l) - delta(i, k)*delta(j, l) - delta(i, l)*delta(j, k)
    del_fact1 = 4*delta(i, k)*delta(j, l) - delta(i, j)*delta(k, l) - delta(i, l)*delta(j, k)
    del_fact2 = 4*delta(i, l)*delta(j, k) - delta(i, k)*delta(j, l) - delta(i, j)*delta(k, l)
    
    accumulator = 0
    called_indices = []
    for J in coordinates:
        for K in coordinates:
            accumulator += np.add(np.add(del_fact0*gamma[J][J][K][K], del_fact1*gamma[J][K][J][K]),  del_fact2*gamma[J][K][K][J])
            called_indices.append((J,J,K,K))
            called_indices.append((J,K,J,K))
            called_indices.append((J,K,K,J))
    accumulator /= 30 
    
    called_indices_no_rep = []
    for indices in called_indices:
        if not indices in called_indices_no_rep:
            called_indices_no_rep.append(indices)
    
    
    return accumulator, called_indices_no_rep

                          
def alpha_alpha_contribution_to_chi3(alpha1, alpha2, i, j, k, l):
    """
    see Kwak 2015 Rigorous theory of molecular orientational nonlinear optics for details
    quote: "the <alpha alpha> class is the contribution from the anisotropic reorientation of the induced dipole moments"
    """
    coordinates = [0, 1, 2]
    
    #Calculate pure electronic case
    del_fact0 = 4*delta(i, j)*delta(k, l) - delta(i, k)*delta(j, l) - delta(i, l)*delta(j, k)
    del_fact1 = 4*delta(i, k)*delta(j, l) - delta(i, j)*delta(k, l) - delta(i, l)*delta(j, k)
    del_fact2 = 4*delta(i, l)*delta(j, k) - delta(i, k)*delta(j, l) - delta(i, j)*delta(k, l)
    
    accumulator = 0
    called_indices = []
    for J in coordinates:
        for K in coordinates:
            accumulator += np.add(np.add(del_fact0*alpha1[J][J]*alpha2[K][K],
                                         del_fact1*alpha1[J][K]*alpha2[J][K]),
                                         del_fact2*alpha1[J][K]*alpha1[K][J])/30
            accumulator -= alpha1[J][J]*alpha2[K][K]*delta(i, j)*delta(k, l)/9
            called_indices.append((J,J,K,K))
            called_indices.append((J,K,J,K))
            called_indices.append((J,K,K,J))
    
    called_indices_no_rep = []
    for indices in called_indices:
        if not indices in called_indices_no_rep:
            called_indices_no_rep.append(indices)
    
    
    return accumulator, called_indices_no_rep

def make_rot_ave_gamma_tensor(gamma):
    coordinates = [0, 1, 2]
    gamma_ave = [[[[None for l in coordinates] for k in coordinates] for j in coordinates] for i in coordinates]
    for i in coordinates:
        for j in coordinates:
            for k in coordinates:
                for l in coordinates:
                    gamma_ave[i][j][k][l], called_indices_no_rep = rot_ave_gamma(gamma, i, j, k, l)
    return np.array(gamma_ave)

def make_uncorrected_chi3_tensor(gamma, alpha1, alpha2):
    gamma_ave = [[[[None for l in coordinates] for k in coordinates] for j in coordinates] for i in coordinates]
    for i in coordinates:
        for j in coordinates:
            for k in coordinates:
                for l in coordinates:
                    gamma_temp, called_indices_no_rep = rot_ave_gamma(gamma, i, j, k, l)
                    alpha_alpha_temp, called_indices_no_rep = alpha_alpha_contribution_to_chi3(alpha1, alpha2, i, j, k, l)
                    gamma_ave[i][j][k][l] = gamma_temp + alpha_alpha_temp
    return np.array(gamma_ave)

def Lorentz_Lorenz_local_field_correction_NB_and_density(gamma, numb_density, sq_refractive_index, lambda_out, lambda_1, lambda_2, lambda_3):
    """
    gamma is a 4d numpy array
    all lambdas in microns
    """
    chi3_conversion_factor = constants.physical_constants["atomic unit of 2nd hyperpolarizability"][0]/constants.epsilon_0
    
    e_0, e_1, e_2, e_3 = sq_refractive_index(lambda_out), sq_refractive_index(lambda_1), sq_refractive_index(lambda_2), sq_refractive_index(lambda_3) 
    correction = ((e_0 + 2)/3)*((e_1 + 2)/3)*((e_2 + 2)/3)*((e_3 + 2)/3)
    print(correction)
    return np.multiply(correction, gamma)*chi3_conversion_factor*numb_density*6

def display_chi3_elements(chi3_symbols, chi3_values):
    """
    when formating with latex
    {0} = {1} \\text{{ }} \\text{{m}}^2 / \\text{{V}}^2
    """
    coordinates = [0, 1, 2]
    for i in coordinates:
        for j in coordinates:
            for k in coordinates:
                for l in coordinates:
                    print("{0} = {1} m^2 / V^2".format(sp.latex(chi3_symbols[i][j][k][l]), chi3_values[i][j][k][l]))

def initialize_3D_4th_rank_tensor(gamma_tuples):
    coordinates = [0, 1, 2]
    gamma = [[[[None for l in coordinates] for k in coordinates] for j in coordinates] for i in coordinates]
    for indices, gamma_val in gamma_tuples:
        gamma[indices[0]][indices[1]][indices[2]][indices[3]] = gamma_val
    return gamma

def initialize_3D_alpha_matrix(alpha_tuples):
    coordinates = [0, 1, 2]
    alpha = [[None for j in coordinates] for i in coordinates]
    for indices, alpha_val in alpha_tuples:
        alpha[indices[0]][indices[1]] = alpha_val
        alpha[indices[1]][indices[0]] = alpha_val
    return alpha

def compute_and_display_chi3_from_raw_gamma(gamma_tuples, numb_density, sq_refractive_index, lambda_out, lambda_1, lambda_2, lambda_3):
    """
    Parameters:
        gamma_tuples: list of tuples
            gamma value in atomic units
            [(axis0, axis1, axis2, axis3), gamma_value),...]
        lambda_out: float
            output frequency in units of microns
        lambda_1: float
            first input frequency in units of microns
        lambda_2: float
            second input frequency in units of microns
        lambda_3: float
            third input frequency in units of microns
        numb_density: float
            density of target molecule in #/m**3
    """

    chi3_elms = sp.symbols("chi^(3)_x:zx:zx:zx:z")
    coordinates = [0, 1, 2]

    # initializing chi3 symbolic tensor
    chi3_sym = [[[[chi3_elms[27*i + 9*j + 3*k + l] for l in coordinates] for k in coordinates] for j in coordinates] for i in coordinates] 

    #initialize gamma NB
    gamma = initialize_3D_4th_rank_tensor(gamma_tuples)

    #initialize gamma NB rot ave
    gamma_rot_ave = make_rot_ave_gamma_tensor(gamma)

    #calculate coorections and change in units
    chi3_rot_ave = Lorentz_Lorenz_local_field_correction_NB_and_density(gamma_rot_ave, numb_density, sq_refractive_index, lambda_out, lambda_1, lambda_2, lambda_3)

    # display results of chi3
    #display_chi3_elements(chi3_sym, chi3_rot_ave)
    
    return gamma_rot_ave, chi3_rot_ave, chi3_sym
                          
def compute_and_display_chi3_from_raw_gamma_alpha(gamma_tuples, alpha_tuples, numb_density, sq_refractive_index, lambda_out, lambda_1, lambda_2, lambda_3):
    "assumes all lambdas are the same"
    chi3_elms = sp.symbols("chi^(3)_x:zx:zx:zx:z")
    coordinates = [0, 1, 2]

    # initializing chi3 symbolic tensor
    chi3_sym = [[[[chi3_elms[27*i + 9*j + 3*k + l] for l in coordinates] for k in coordinates] for j in coordinates] for i in coordinates] 

    #initialize gamma
    gamma = initialize_3D_4th_rank_tensor(gamma_tuples)

    #initialize alpha
    alpha = initialize_3D_alpha_matrix(alpha_tuples)

    #make uncorrected chi3 tensor
    chi3_uncorrected = make_uncorrected_chi3_tensor(gamma, alpha, alpha)

    #calculate coorections and change in units
    chi3_corrected= Lorentz_Lorenz_local_field_correction_NB_and_density(chi3_uncorrected, numb_density, sq_refractive_index, lambda_out, lambda_1, lambda_2, lambda_3)

    # display results of chi3
    #display_chi3_elements(chi3_sym, chi3_corrected)
    
    return chi3_uncorrected, chi3_corrected, chi3_sym
