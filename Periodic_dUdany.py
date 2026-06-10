# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 21:31:57 2026

@author: DELL
"""

import numpy as np
#import sympy as sp
import scipy.integrate as spint
from functools import reduce
from scipy.linalg import expm
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from collections import Counter
import pandas as pd
from pathlib import Path

#%% code up winding number due to bulk hamiltonian 

def J_tilte_a(t, t_a, pulse_length, l, T, pulse_amp): # t_a define the START POINT in each period T for pulse a
    t_p = np.mod (t, T) #for periodic changings, but there's problem for Jz needed to be figured out

    t_start_a = l * T + t_a
    
    t_end_a = t_start_a + pulse_length
    
    if t_end_a < T:
    
#    if t_start_a < t < t_end:        
#        return  pulse_amp   
#    else:
#        return 0 -- this if else is not suitable when t is an array

        return np.where((t_p >= t_start_a) & (t_p < t_end_a), pulse_amp, 0)
    #np.where(condition, x, y) if condition is true : take x(pulse amp here) ; 
    #                          if condition is false: take y (0)3.+1
    if t_end_a > T:
        
        t_end_a_new = t_end_a - T
        
        return np.where((t_p >= t_start_a) & (t_p < T) | (t_p >= 0) & (t_p < t_end_a_new), pulse_amp, 0)
    
#parameter adjustion:
#(here always take T = 1, l = 0, set the numerical platform of time, also unify the y axis as epsilon in fig_1)
T = 1
l = 0
#(as said in the paragraph above (16), t_x = 0, t_y = T/3, t_z = 2*T/3)
t_x = 0
t_y = T/3
t_z = T * (2/3)
    
#delta_t = T/2
delta_t = T/2

J_a = 2* ((0.9 * np.pi ))/ ((T * 0.5) * 4)
#J_a = 2* ((0.9 * np.pi ))/ ((T * 0.5) * 4)  #  2* ((0.9 * np.pi ))/ ((T * 0.5) * 4)-- this is 2x the value in paper 
#J_a =(1* np.pi)/3# try new J_a
                                       # Q: Can I vary this J_a to perform perturb on Hamiltonian?
                                       # Ans: I don't think so. I think this might work: 1.change the pulse amplitude by 1%
                                       #                                                 2.shift the timing by 0.01T
                                       #                                                 3.add a tiny disorder term


#%% (here change H_test to real H, which is parameterised by k:
#(firstly try to use the [H] build a list of 6 H; 
#(another method is that to firstly write H1/U1 to H6/U6, then U_full = U1U2U3...U6)


# QUESTION: should I define a function: array of Hamiltonians? known H_i depend on Jx, Jy, Jz, k
# -- firstly we try only do as quantities
# here t_x,y,z are setted, delta_t are setted, J_x,y,z are setted 
def J_tilte_x(t):
    return J_tilte_a(t, t_x, delta_t, l, T, J_a)
def J_tilte_y(t):
    return J_tilte_a(t, t_y, delta_t, l, T, J_a)
def J_tilte_z(t):
    return J_tilte_a(t, t_z, delta_t, l, T, J_a) 
# didn't write a loop here because J_a and J_start are also variables, can be different/change 
# we defined fixed l = 1, T_length = T = 1 
#%%
# try plot: t is from 0 to 1, J_tilte_x, y, z 

#print (J_tilte_x(0))
#print (J_tilte_x(0.5))

#print (J_tilte_y(1/3))
#print (J_tilte_x(5/6))
# here is clear that each pulse is defined as start point(of driving) = Ja, end point (of driving) = 0
#%%
t_test = np.linspace(-2, 4, 100) 
# function of J_a need to be adjusted as in that t can be out of range of one period(0,T)

plt.plot(t_test, J_tilte_x(t_test), label='J_tilte_x')
plt.plot(t_test, J_tilte_y(t_test), label='J_tilte_y')
plt.plot(t_test, J_tilte_z(t_test), label='J_tilte_z')

plt.xlabel('t')
plt.ylabel('J_tilte_a')
plt.title('Plot of J_tilta_a(pulse of hopping amp) vs t')
plt.legend()
plt.grid(True)
plt.show()

#%% some test try figure properties of J func
#t_specific_time_point_test = np.pi - 3

#call_func_J_tilte_x = J_tilte_z(t_specific_time_point_test)

#print (call_func_J_tilte_x, ' see if J_tilte_a = x return a scalar multiply the Pauli matrices')
#plot succeed

#%%

# **a_0 is lattice spacing**, which is the quantity ensures translational invariance(in x direction):
# H(r+ a_0 *x_hat, r' + a_0 *x_hat, t) = H(r, r', t) -- here H is Hamiltonian and x_hat is unit vector in x direction

# Q:1.  what is r, r'(both vectors) specifically here? do they refer to interactions between r and r'? 
# Q:2.  does (r-r' )'s modulus<l shows some kind of interaction range between two sites/positions? 
#(just like the interaction in SSH model, but in SSH model H11 doesn't represent real position(as in meter or...) - )

#(- it represent onsite interaction: hopping amplitude between site 1 and site 1 )
# a_0 is also called 'lattice_length' in the subsequent codes, (coz a_0 was used directly in parameter in functions?

#lattice_length = 1
a_0 = 1
#stage_length = T/6 #(numerical calculations)
# due to equation (28)
#stage_number = 6 #(just N afterwards)
#Pauli matrices:
pauli_x = np.array([[0.0, 1.0], [1.0, 0.0]])
pauli_y = np.array([[0.0, -1j], [1j, 0.0]])
pauli_z = np.array([[1.0, 0.0], [0.0, -1.0]])

# Pulse setting Q: where is J_0? Please find it later
    
#%%
#-----K BASIS HAMILTONIAN MATRIX H 2X2-----
#Here is the 2by2 Hamiltonian that describes the [sublattices coupling] in momentum space
#            -- we will refer this as '2by2 Hamiltonian' /'2by2 H' in the future

#delta_t is the actual T/2, pulse length of our case 
def H_matrix_2by2(k_x, k_y, t, a_0): 
    
    def J_tilte_x(t):
        return J_tilte_a(t, t_x, delta_t, l, T, J_a)
    def J_tilte_y(t):
        return J_tilte_a(t, t_y, delta_t, l, T, J_a)
    def J_tilte_z(t):
        return J_tilte_a(t, t_z, delta_t, l, T, J_a) 
    
    H = (J_tilte_x(t) + J_tilte_y(t) * np.cos (k_x * a_0) + J_tilte_z(t) * np.cos(k_y * a_0)) * pauli_y + (J_tilte_y(t) * np.sin (k_x * a_0) + J_tilte_z(t) * np.sin(k_y * a_0)) * pauli_x
    #(28)
    
    return H

# generate a list of t which takes 6 samples of static H(1 to 6)

#t = T/6 * i + T/12 (T/12 is for if there is any confusion in the boundary points of the six stages)


# run a small test for the future 2N 2N matrix
#H_test_for_entrices = H_matrix_2by2(np.pi/2, np.pi/2, 5/12, 1)
#a_11 = H_test_for_entrices[0][0]
#a_12 = H_test_for_entrices[0][1]
#(1.413716694115407-1.413716694115407j) (1.413716694115407-1.413716694115407j)

#%% Guess: bulk-defect as in bulk hamiltonian's winding number decide what will happen when we introduce defect
#hence we don't need to calculate winding number through the one with time vortex; we can try: without time vortex
# write H_eff -- H_eff is calculated from U(k, T)

def U_full_discretisation(k_x, k_y, num_time_stages): # which means along lhs to rhs, from (1) earlier/smaller t to (6) later/larger t

    U_list_N_matrices = []
    t_desicretisation= T /num_time_stages
#    time_stage = num_time_stages
    
    for i in range (num_time_stages):
        
    
        t_i = T - t_desicretisation * (i + 1) + t_desicretisation/2 # here shift sample t_i to midpoint 
                                                                    # i = 0 (leftest in muti) t_i = T - dt + dt/2
                                                                    # i = N -1 t_i = T - (T/N) * N + T/2N
                                                                    # this will be safe for any case
        H_i = H_matrix_2by2(k_x, k_y, t_i, a_0)
#expm(-1j* H_i_test* delta_t_test) 
        
        U_i = expm(-1j* H_i* (t_desicretisation)) #!!!!! # here using exponential form construct each of U_1 to U_6
        
        U_list_N_matrices.append(U_i)
        
#        print(U_list_N_matrices)
        
    U_full = reduce(np.matmul, U_list_N_matrices) # here adjoint the full time evolution 
        # here U_list_N_matrices = [e^{-iH(T-dt)}]
    return U_full   
# this will be for U(k, t)    

#%% this is for not U(k, T)/U(r, T), but U(k/r, t)
def U_t_discretization(k_x, k_y, num_time_stages, t): # which means along lhs to rhs, from (1) earlier/smaller t to (6) later/larger t

    U_list_N_matrices = []
    t_discretization= t /num_time_stages
#    time_stage = num_time_stages
    
    for i in range (num_time_stages):
        
    
        t_i = t - t_discretization * (i + 1)  + t_discretization/2   #t_discretization = delta_t here, **take midpoint  
        # here we shifted to midpoint 
        # note: here t_discretization << T/6
        # because if t is any value, and no classification of t range is done ..|.....|.....|.....| 
        # t/ N_dis doesn't need to not across the '|' and blur the change from Hi to Hi+1 at the boundaries 
        # double check here: i = 0 (leftest), t_i = t - delta_t
     #                       i = num - 1(rightest) , t_i = t - num * dis = t - t = 0 
        H_i = H_matrix_2by2(k_x, k_y, t_i, a_0)
#expm(-1j* H_i_test* delta_t_test) 
        
        U_i = expm(-1j* H_i* (t_discretization)) #!!!!! # here using exponential form construct each of U_1 to U_6
        
        U_list_N_matrices.append(U_i)
        
#        print(U_list_N_matrices)
        
    U_full = reduce(np.matmul, U_list_N_matrices) # here adjoint the full time evolution 
        # here U_list_N_matrices = [e^{-iH(T-dt)}]
    return U_full   
# this will be for U(k, t)
#%%
#test_reduce()
#test_x = np.array([[1j +1, 3j -2], [1, 3]])
#test_y = np.array([[1j + 0.5, 0.2j + 5], [5, 1]])

#list_test = [test_x, test_y]
#list fruit = [apple, pear, peach]
#U_test = reduce(np.matmul, list_test) 
#print(U_test)
# all verified 
#%% here is a simple version can be applied to U_full_k as in 6 constant stages

def U_full_descending(k_x, k_y, num_time_stages): # U_full =  U6 U5 U4 U3 U2 U1 as in U1 is from 0 to T/6

# difference from discretisation func: H(t) sample take t in the midpoint of H const stage 
    U_list_N_matrices = []
    t_desicretisation= T /num_time_stages
#    time_stage = num_time_stages
    
    for i in range (num_time_stages):
        
    
        t_i = T - t_desicretisation * (i + 1/2) # from last one to first one
    
        H_i = H_matrix_2by2(k_x, k_y, t_i, a_0)
#expm(-1j* H_i_test* delta_t_test)
        
        U_i = expm(-1j* H_i* (t_desicretisation)) # H_0 will be H (t = T-(T/6)*1/2), with is 6th stage
        
        U_list_N_matrices.append(U_i)
        
    U_full = reduce(np.matmul, U_list_N_matrices) # here adjoint the full time evolution 
        # here U_list_N_matrices = [e^{-iH(T-dt)}]
    return U_full 

#%%
#TEST : if U_full is unitary? 
#U_full_6TS_wrong =  U_full_discretisation(np.pi/2, np.pi/5, 6) # these kx, ky can be changed 
# shouldn't risk to use discretisation func to calculate the 6 const H partern, for that the start/end point of a pulse 
# might be defined as 0 or Ja, uncertain

#adjoint of U_full_6TS(U_adgga)
#U_full_6TS_wrong_adj = U_full_6TS_wrong.conj().T
#print("\nU† =", U_full_6TS_wrong_adj)
# check if U is really 'unitary'
#check_unitariness = U_full_6TS_wrong @ U_full_6TS_wrong_adj
#print('unitary?', check_unitariness) # -- we should see this is always identity(U is indeed Unitary)
#ANS： Yes. U_full (U(k, T)) is unitary
#But:  It doesn't eqaul to identity at all; so due to RLBL article III-A, it needs to be modified 
#Q: is U diagonalised? If so it would be much more better! Since we don't need to include 'V' as in U_dia = V U V^-1 subsequently
#%% the logarithm is extracting the eigenphases of U when U is unitary 
# H_eff = template code -- 
def heff_from_U(U, T, epsilon, tol=1e-10):
    """
    Compute H_eff = (i/T) log_epsilon(U)
    with the logarithm branch cut along exp(-i * epsilon * T).

    Parameters
    ----------
    U : (N,N) complex ndarray
        Floquet/unitary operator.
    T : float
        Period.
    epsilon : float
        Quasienergy value defining the branch cut.
    tol : float
        Numerical tolerance for unitarity check.

    Returns
    -------
    Heff : (N,N) complex ndarray
        Effective Hamiltonian.
    evals : (N,) complex ndarray
        Eigenvalues of U.
    quasienergies : (N,) float ndarray
        Branch-chosen quasienergies used in Heff.
    """
    U = np.asarray(U, dtype=complex)
    N = U.shape[0]

    # Optional: check unitarity
    I = np.eye(N, dtype=complex)
    if not np.allclose(U.conj().T @ U, I, atol=tol):
        print("Warning: U is not exactly unitary within tolerance.")

    # Diagonalize U
    evals, V = np.linalg.eig(U)

    # For a unitary matrix, eigenvalues should lie on the unit circle:
    # evals_n = exp(i * alpha_n), where alpha_n is an angle.
    alpha = np.angle(evals)   # principal angles in (-pi, pi]

    # Branch cut direction is exp(-i * epsilon * T)
    cut = -epsilon * T

    # Put alpha on the interval (cut - 2*pi, cut]
    alpha_branch = ((alpha - cut) % (2 * np.pi)) + cut
    alpha_branch = np.where(alpha_branch > cut, alpha_branch - 2 * np.pi, alpha_branch)

    # Since log(lambda_n) = i * alpha_branch_n for |lambda_n|=1,
    # Heff eigenvalues are:
    quasienergies = -alpha_branch / T

    # Reconstruct Heff
    V_inv = np.linalg.inv(V)
    Heff = V @ np.diag(quasienergies) @ V_inv

    # Clean up tiny non-Hermitian numerical noise
    Heff = 0.5 * (Heff + Heff.conj().T)

    return Heff, evals, quasienergies

#%% write my own Heff code

# find V as similarity trans matrix. U_diagonalised = V U V^-1 = (diagonal: eigenval = e^(-j Phi_1)... e^(-j n))

#U_full_bulk is the 2x2 full time evolution opt difined by different kx, ky, but not t as in they are U_full = U (k_vec, T)

# later maybe can embed the U_full_bulk(T) = [e ^ (-j *H(k,t_n_last) *delta_t)] ... [e ^ (-j *H(k,t_n_first) *delta_t)] 
# where H2x2 calculated by function H_matrix_2by2, with input k_x, k_y, t_n, a_0(as fixed lattice length)
# from INPUT: kx, ky, t 

# Q1: what is the upper and lower limitation for the winding number integeral(equation (4)) when U is general case?
# Q2: does general case U_epsilon(k, t) we setted, use the same definition as U(k, t) as in t_vor article, eqn 21 (U(...t) does mean to multiply all H up to t)
#     - Ans: the first half in one period (T) have U(k, 2t) used same; hence should use func H(k,t) >> times up to U(t)
#            double check the maths behind why we can (ex:U(T) = U6U5...U1) break U(t)= U(t)...U(1)
#     - Ans: the second half in one period (T) = V_epsilon, which is from U(T)'s corresponding effective Ham Heff 


#LATER ADAPT U_full to: when calculating Heff, INPUT needed now: U(T); but we can also break U(full) into that we define H(k, t) and t descritisation(summed over)
# INPUT will be: kx, ky, num of discritisation(for T/n = delta t, H(t) in expoential of U_n...U_1 is H( T/n * i)); subQ: can we introduce time vortex in bulk k space Hamiltonian? 
#                                                                                                                        If we cannot, why? 


#%% Rewrite the function above: 
#   INPUT: U_full_bulk = U(k, T) epislon -- which is used to define the branchcut

def Heff_from_Ufull (U_full_bulk, epsilon): # U_full_bulk = U(k, T)
    
    eigvals_Ufull, eigvecs_Ufull = np.linalg.eig(U_full_bulk)  #Each eigval_alpha corresponding to each eigenvector |phi_alpha>
    
    U_dia_full = np.diag(eigvals_Ufull) # diagonalised U, now each element on the diagonal is lamda_alpha = eigval_alpha corresponding to eigenvector |phi_alpha>
    
    V = eigvecs_Ufull # need to be checked that if they are natrually returned in columns  # ***VERIFICATION 1
    
    PHI_dia = np.angle(U_dia_full) #np.angle works elementwise on each entry and returns the phase angle of each complex number, in the range (-pi, pi]
    #take angle as an array maintian the order of alpha(s) perfectly
    #**VERIFICATION2**: U_diag = np.exp(1j*PHI_diag) instead of **elementwise** definition of phi_alpha = - PHI_alpha

    #Q: expm(1j*M) is matrix exponential -- how does it work diffrently from elementwise exp? 
    PHI_alpha_array = np.diagonal(PHI_dia) # branch range: (-pi, pi]



    PHI_modified_list = []
    
    # possible modified branchrange: 
        #epsilonT range: [-pi, pi] PHI_modified_min > -epsilonT_min - 2pi = -3pi; PHI_modified_max <= -espilonT_max = -(-pi) = pi 
    # hence (rough) PHI_modified range (with above information): (-3pi, pi]
    # we only test PHI_alpha minus 2npi, because only the smallest PHI, (-pi), + 2pi will be in the range of PHI_modified;
    # That need epsilon T to be -pi, hence - epsilon * T = -(-pi) = pi; - epsilon * T - 2 * np.pi = -pi, PHI_alpha = -pi is out of range--> modified to pi
    # but this is a special case, 
    
    for PHI_alpha in PHI_alpha_array:
        
        if  - epsilon * T - 2 * np.pi < PHI_alpha <= - epsilon * T:   # here used (] for PHI, as we take -epsilonT but not  - epsilonT - 2* pi
            
            PHI_modified = PHI_alpha
        
        elif  - epsilon * T - 2 * np.pi <  PHI_alpha - 2 * np.pi <=  - epsilon * T:
            
            PHI_modified = PHI_alpha - 2 * np.pi
        
        elif - epsilon * T - 2 * np.pi <  PHI_alpha - 4 * np.pi <=  - epsilon * T:
            
            PHI_modified = PHI_alpha - 4 * np.pi
        
        elif - epsilon * T - 2 * np.pi <  PHI_alpha + 2 * np.pi <=  - epsilon * T:
            
            PHI_modified  =  PHI_alpha + 2 * np.pi
            print(PHI_alpha, PHI_modified) # this step is to check if the guess, only case here should be -pi + 2pi = pi, is correct
            
        else:
            raise ValueError(f"No condition matched for PHI_alpha = {PHI_alpha}")

        
        # also remember: in one branch each PHI_alpha_modified only has one value choice
        PHI_modified_list.append(PHI_modified)
        
    phi_alpha_array = - np.array(PHI_modified_list)
    
    H_eff_dia = (1/ T) * np.diag(phi_alpha_array)

    H_eff = V @ H_eff_dia @ np.linalg.inv(V)

    return H_eff
#VERIFICATION 2: within the code above 
#VERIFICATION 3: pen and paper: what if we go from phi_diag = - np.angle (U_dia_full)? Done
#VERIFICATION 4: input a real matrix Phi to calculate Heff, see if there's any error.Also compare both function Heff_from_Ufull and Heff_from_Ufull 2 output
#VERIFICATION 5: test if V is as what expected

#%%
def VERIFICATION2(U_full_bulk):
    
    eigvals_Ufull, eigvecs_Ufull = np.linalg.eig(U_full_bulk)  #Each eigval_alpha corresponding to each eigenvector |phi_alpha>
    
    U_dia_full = np.diag(eigvals_Ufull) # diagonalised U, now each element on the diagonal is lamda_alpha = eigval_alpha corresponding to eigenvector |phi_alpha>

    PHI_diag = np.angle(U_dia_full)
    phi_diag = - np.angle(U_dia_full)

    
    if np.allclose(expm(1j * PHI_diag), U_dia_full):
        
        if np.allclose(expm(-1j * phi_diag), U_dia_full):
            
         print ('VERIFIED')
        else:
            
            raise ValueError('error TEST2')
    
    else:
        raise ValueError('error TEST1')
        
    return 

#F = VERIFICATION2(U_full_descending((2*np.pi)/8, (2*np.pi)/8, 6))
# VERIFICATION2: verified 

#%%

def Heff_from_Ufull_2 (U_full_bulk, epsilon): # Heff_from_Ufull_2 has no PHI, only phi_alpha 
    
    eigvals_Ufull, eigvecs_Ufull = np.linalg.eig(U_full_bulk)  #Each eigval_alpha corresponding to each eigenvector |phi_alpha>
    
    U_dia_full = np.diag(eigvals_Ufull) # diagonalised U, now each element on the diagonal is lamda_alpha = eigval_alpha corresponding to eigenvector |phi_alpha>
    
    V = eigvecs_Ufull # need to be checked that if they are natrually returned in columns  # ***VERIFICATION 1
    
    phi_dia = - np.angle(U_dia_full) #[-pi, pi)
    #take angle as an array maintian the order of alpha(s) perfectly

    phi_alpha_array = np.diagonal(phi_dia) # branch range: (-pi, pi]
    
    phi_modified_list = []
    
    # possible modified branchrange: 
        #epsilonT range: [-pi, pi] PHI_modified_min > -epsilonT_min - 2pi = -3pi; PHI_modified_max <= -espilonT_max = -(-pi) = pi 
    # hence (rough) PHI_modified range (with above information): (-3pi, pi]
    # we only test PHI_alpha minus 2npi, because only the smallest PHI, (-pi), + 2pi will be in the range of PHI_modified;
    # That need epsilon T to be -pi, hence - epsilon * T = -(-pi) = pi; - epsilon * T - 2 * np.pi = -pi, PHI_alpha = -pi is out of range--> modified to pi
    # but this is a special case, 
    
    for phi_alpha in phi_alpha_array:
        
        if  epsilon * T <= phi_alpha <  epsilon * T + 2 * np.pi :   # here used (] for PHI, as we take -epsilonT but not  - epsilonT - 2* pi
            
            phi_modified = phi_alpha
        
        elif  epsilon * T <= phi_alpha + 2 * np.pi <  epsilon * T + 2 * np.pi:
            
            phi_modified = phi_alpha + 2 * np.pi
        
        elif  epsilon * T <= phi_alpha + 4 * np.pi <  epsilon * T + 2 * np.pi:
            
            phi_modified = phi_alpha + 4 * np.pi
        
        elif epsilon * T <= phi_alpha - 2 * np.pi <  epsilon * T + 2 * np.pi:
            
            phi_modified = phi_alpha - 2 * np.pi
            print(phi_alpha, phi_modified) # this step is to check if the guess, only case here should be -pi + 2pi = pi, is correct
            
        else:
            raise ValueError(f"No condition matched for phi_alpha = {phi_alpha}")

        
        # also remember: in one branch each PHI_alpha_modified only has one value choice
        phi_modified_list.append(phi_modified)

    
    H_eff_dia = (1/ T) * np.diag(np.array(phi_modified_list))

    H_eff = V @ H_eff_dia @ np.linalg.inv(V) # A *B is elementwise multiplication but not matrix multiplication 

    return H_eff

#%%

# Which U equation is recalled? U_full_descending, which is only 6 stages( num_of_stages = 6, fixed, hence t_discretisation = T/6 which is large)
#                                                                          hence here the H_i is taken as the mid-point of each t_discre interval
#                               U_discretisation input with large N, very small t_discretisation; hence H_i is taken as lower edge point (0), delta_t, ...
def Heff_from_kxky (k_x, k_y, num_time_stages, epsilon):
    U_kxky = U_full_descending(k_x, k_y, num_time_stages)   # U(k, T)
    Heff = Heff_from_Ufull (U_kxky, epsilon)
    return Heff

def Heff_from_kxky_2(k_x, k_y, num_time_stages, epsilon):
    U_kxky = U_full_descending(k_x, k_y, num_time_stages)   # U(k, T)
    Heff = Heff_from_Ufull_2 (U_kxky, epsilon)
    return Heff

# here for bulk as in infinite amount of sites, or infinite length of system, no edge, N --> inf; k_x, k_y spacing -->0

def VERIFICATION_3(k_x, k_y, epsilon):
    Heff_1 = Heff_from_kxky (k_x, k_y, 6, epsilon)
    Heff_2 = Heff_from_kxky_2 (k_x, k_y, 6, epsilon)
    
    if np.allclose(Heff_1, Heff_2):

         print ('VERIFIED')
    else:
            
        raise ValueError('error TEST')
        
    return 

#%% can be commented out
#Heff_test = Heff_from_kxky ((2*np.pi)/8, (2*np.pi)/7, 6, np.pi/2)

#U_full_test = U_full_descending((2*np.pi)/8, (2*np.pi)/7, 6)

#Heff_test_exp = expm(-1j * Heff_test * T)
#TEST1 = np.allclose(Heff_test_exp, U_full_test)
#print(TEST1)

#%% can be commented out
#Heff_test_2 = Heff_from_kxky_2 ((2*np.pi)/8, (2*np.pi)/7, 6, np.pi/2)

#Heff_test_exp_2 = expm(-1j * Heff_test_2 * T)
#TEST2 = np.allclose(Heff_test_exp_2, U_full_test)
#print(TEST2)
#%% can be commented out
#Heff_verify = VERIFICATION_3((2*np.pi)/8, (2*np.pi)/7, 0)
    
#%%
#eigvals_Ufull, eigvecs_Ufull = np.linalg.eig(U_full_test)  #Each eigval_alpha corresponding to each eigenvector |phi_alpha>
    
#U_dia_full = np.diag(eigvals_Ufull) # diagonalised U, now each element on the diagonal is lamda_alpha = eigval_alpha corresponding to eigenvector |phi_alpha>
    
#V = eigvecs_Ufull

#U_return = V @ U_dia_full @ np.linalg.inv(V) #wrong: U_diagonalised = V U V^-1; true:  A=VDV−1
#TEST3 = np.allclose(U_return, U_full_test)
#print(TEST3)

#%%

#equation 3 as in the definition of U(k, t) Done; U(k, t) = U_k (t, 0)

#I equa 8 finish   

#II integral   

#%% Equation 8

# the U_e function below can only calculate U(t) with t within [0, T]
# only with H(2x2) for kx, ky, without time vor
# -- otherwise : H_i = H2NyNx_xyOBC_time_vor_core (N_x, N_y, t_i, a_0, x_0, y_0) instead of H_matrix_2by2


def U_epsilon_H2by2_allt (kx, ky, t, N_dis, epsilon): # N_dis stand for time discretisation we choose to calculate U(k,t) here 
                                 #only_bulk_U_involved will be =True or = False
                                 #Physically 'only_bulk_U_involved' means no time vortices 

    U_full_for_Heff = U_full_discretisation(kx, ky, num_time_stages = 6) # for full period U, 
                                                                         # H2x2, can be divided into 6 stages piecewise H 
                                                                         # each H takes the midpoint of that stage 
    #**if any change on H, H isn't valid for time discretization only 6 (piecewise) anymore, above should be changed 
    
    
    if 0<= t <= T/2:
        
        t_epsilon = 2 * t
    # alternative way record(little notebook p79): we can just not call any U function before, but classify t_epsilon into 6 stages, for each stage we have a certain parttern of 
        U_epsilon =  U_t_discretization(kx, ky, N_dis, t_epsilon) #here the discretization of time is done only in the first half 
                                       # here was set as 6 in the could we send to Mimi intially, wrong!
                    #U_t_discretization(k_x, k_y, num_time_stages, t) - these are meaning of variables 
       #U_e (k, t) = U(k, 2t) = U_latest......U_earliest which should be calculated in numerical(infinitesimal) pieces
        # for efficiency: 
        
    elif T/2 <= t <= T:  # here this must be elif, otherwise all t belongs to the first 'branch' 
                         # will finally go to else, which raise error interrupt the function
        
        t_epsilon = 2 *T - 2 * t # 2nd def
        
        H_eff_epsilon = Heff_from_Ufull_2(U_full_for_Heff, epsilon) 
        #Here this H_eff_epsilon has specific k_vec as U_full_for_Heff
        # here used the second definition (with phi, as - whatever on the phase part of a complex num, but not PHI)
        U_epsilon =  expm(-1j* H_eff_epsilon * t_epsilon ) # H_eff_epsilon refer to above, t_epsilon as 2nd def
        
    else:
        raise ValueError(f"No condition matched for U(k_vec, t) with time  = {t}")
    return U_epsilon  


#%%
def H_epsilon_H2by2_first (kx, ky, t, epsilon): 
    # this doesn't produce U, but only produce H(t_[some index])
    # e^H(t_[some index])*delta_t will be followed, time discretization will depend on index but not how t is input here
    #Note: here we don't sample t 

        
    t_epsilon = 2 * t
    H_22_kxkyt =  H_matrix_2by2(kx, ky, t, a_0)
       
    return H_22_kxkyt
    
def U_epsilon_H2by2_second (kx, ky, t, epsilon):
       
    U_full_for_Heff = U_full_discretisation(kx, ky, num_time_stages = 6)
    t_epsilon = 2 *T - 2 * t # 2nd def
        
    H_eff_epsilon = Heff_from_Ufull_2(U_full_for_Heff, epsilon) 
    U_epsilon =  expm(-1j* H_eff_epsilon * t_epsilon ) # H_eff_epsilon refer to above, t_epsilon as 2nd def

    return U_epsilon  

def periodic_derivative(arr, spacing, axis):
    """
    Central finite difference with periodic boundary conditions:
        (f[i+1] - f[i-1]) / (2 spacing)
    """
    return (np.roll(arr, -1, axis=axis) - np.roll(arr, 1, axis=axis)) / (2 * spacing)


def winding_number_enumerate_fixed(n_kx=40, n_ky=40, n_t=80, epsilon= 0 ):
    """
    Compute Rudner winding number W[U_epsilon].

    axis 0 -> kx
    axis 1 -> ky
    axis 2 -> t

    This version:
      1. builds one full U_epsilon(kx, ky, t) array,
      2. uses endpoint=False for all periodic variables,
      3. uses physical derivatives d/dkx, d/dky, d/dt,
      4. multiplies by dkx*dky*dt at the end.
    """

    if n_t % 2 != 0:
        raise ValueError("Please use even n_t, e.g. n_t=40, 80, 100.")

    CONST = 8 * np.pi**2

    # Periodic variables: do not include both 0 and 2pi / T.
    kx_list = np.linspace(0, 2*np.pi, n_kx, endpoint=False)
    ky_list = np.linspace(0, 2*np.pi, n_ky, endpoint=False)
    times = np.linspace(0, T, n_t, endpoint=False)

    dkx = kx_list[1] - kx_list[0]
    dky = ky_list[1] - ky_list[0]
    dt = times[1] - times[0]

    n_half = n_t // 2

    # Full storage for U_epsilon(kx, ky, t)
    u = np.zeros((n_kx, n_ky, n_t, 2, 2), dtype=np.complex128)

    for i, kx in enumerate(kx_list):
        print(f"Building U grid: kx index {i+1}/{n_kx}")

        for j, ky in enumerate(ky_list):

            # ============================================================
            # First half:
            # U_epsilon(k,t) = U(k, 2t), 0 <= t < T/2
            # ============================================================
            U_current = np.eye(2, dtype=np.complex128)

            for m in range(n_half):
                t = times[m]

                # Store first. Therefore u[...,0] = identity.
                u[i, j, m, :, :] = U_current

                # Advance from t to t+dt in epsilon-time.
                # Since U_epsilon(t)=U(2t), original physical time step is 2*dt.
                tau_mid = 2 * (t + dt/2)

                H_mid = H_matrix_2by2(kx, ky, tau_mid, a_0)
                U_step = expm(-1j * H_mid * (2*dt))

                U_current = U_step @ U_current

            # Now U_current is the full-period U(k,T) generated consistently
            # with the same discretization used above.
            U_full = U_current

            # Use your template Heff function; it returns (Heff, evals, quasienergies)
            H_eff, _, _ = heff_from_U(U_full, T, epsilon)

            # Optional check:
            if not np.allclose(expm(-1j * H_eff * T), U_full, atol=1e-7):
                print("Warning: exp(-i Heff T) does not match U_full at", i, j)

            # ============================================================
            # Second half:
            # U_epsilon(k,t) = exp[-i H_eff(k) (2T - 2t)],
            # T/2 <= t < T
            # ============================================================
            for m in range(n_half, n_t):
                t = times[m]
                tau = 2*T - 2*t
                u[i, j, m, :, :] = expm(-1j * H_eff * tau)

    # Check periodic closure: U_epsilon(0)=I, U_epsilon(T)=I.
    # T is not included in the grid, but the last point should be close to I
    # as n_t increases.
    print("Max deviation of u(t=0) from I:",
          np.max(np.abs(u[:, :, 0] - np.eye(2))))

    # Inverse
    u_inv = np.linalg.inv(u)

    # Physical derivatives
    dU_dkx = periodic_derivative(u, dkx, axis=0)
    dU_dky = periodic_derivative(u, dky, axis=1)
    dU_dt  = periodic_derivative(u, dt,  axis=2)

    A_kx = np.matmul(u_inv, dU_dkx)
    A_ky = np.matmul(u_inv, dU_dky)
    A_t  = np.matmul(u_inv, dU_dt)

    comm = np.matmul(A_kx, A_ky) - np.matmul(A_ky, A_kx)

    integrand = np.trace(np.matmul(A_t, comm), axis1=-2, axis2=-1)

    W_raw = np.sum(integrand) * dkx * dky * dt
    W = W_raw / CONST

    print("\nRaw integral =", W_raw)
    print("W =", W)
    print("W.real =", W.real)
    print("W.imag =", W.imag)

    return W, u, integrand

W_enum, u_grid, integrand_grid = winding_number_enumerate_fixed(
    n_kx=30,
    n_ky=30,
    n_t=30,
    epsilon= np.pi
)
#%% convergence test
for n in [10, 20, 30, 40]:
    W_enum, _, _ = winding_number_enumerate_fixed(
        n_kx=n,
        n_ky=n,
        n_t=2*n,
        epsilon=np.pi
    )
    print("n =", n, "W =", W_enum)
