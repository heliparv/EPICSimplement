import numpy as np
import pandas as pd
from copy import deepcopy

#get monoculture data
monoculture_abundances = pd.read_excel('Generalized_EPICS/ansari_et_al_average.xlsx', sheet_name=0, header=None).to_numpy()

#number of species
n = len(monoculture_abundances[0])

#get leave-one-out data
loo_abundances = pd.read_excel('Generalized_EPICS/ansari_et_al_average.xlsx',sheet_name=1, header=None).to_numpy()

interactions_matrix = np.zeros((n,n))

#loop calculates interaction terms for one species at a time
for i in range(0,n):
    #make copy of monoculture abundances checking that abundance for target nonzero
    monoc_abun = monoculture_abundances[monoculture_abundances[:,i] != 0.0][:,i]

    #create array that can be concatenated to leave-one-out data
    monoc = np.zeros((len(monoc_abun), n))
    monoc[:, i] = monoc_abun

    #create matrix A for Ax=b to solve for x
    species_matrix = np.r_[monoc, loo_abundances[loo_abundances[:,i] != 0.0]]

    #solves Ax=b for x, creates interaction terms
    interactions_vector = np.linalg.lstsq(species_matrix, -1*np.ones((len(species_matrix),1)), rcond=None)

    interactions_matrix[i, :] = np.reshape(interactions_vector[0], n)

print("Interactions matrix:")
print(interactions_matrix)

#normalizes interaction matrix by self-interactions like done in EPICS paper
normalized_interactions_matrix = deepcopy(interactions_matrix)
for i in range(0, n):
    for j in range(0, n):
        normalized_interactions_matrix[i][j] = normalized_interactions_matrix[i][j]/interactions_matrix[j][j]
print("\nNormalized interactions matrix:")
print(normalized_interactions_matrix)

#Calculates prediction for abundances in n-member community
abundance_n_member_community = np.linalg.lstsq(interactions_matrix, -1*np.ones((n,1)), rcond=None)

print("\nPredicted abundances in an n-member community")
print(abundance_n_member_community[0])