import numpy as np
from scipy.linalg import solve as slv

#Function of EPICS
#Takes number of species, monoculture abundances and abundances in leave-one-out cultures as input
#returns calculated effective pairwise interactions and calculated abundances in n-member community

def epics(monoc, loo):
    abundance_in_monoc = monoc

    n = len(abundance_in_monoc)

    abundance_in_loo = loo

    interactions_matrix = np.zeros((n,n))

    #loop calculates interaction terms for one species at a time
    for i in range(0,n):
        monoc_abun = np.zeros(n)
        monoc_abun[i] = abundance_in_monoc[i]
        
        #create matrix A for Ax=b to solve for x
        species_matrix = np.vstack((monoc_abun, abundance_in_loo[abundance_in_loo[:,i] != 0.0]))

        #solves Ax=b for x, creates interaction terms
        interactions_vector = np.linalg.lstsq(species_matrix, -1*np.ones((len(species_matrix),1)), rcond=None)

        interactions_matrix[i, :] = np.reshape(interactions_vector[0], n)

    #calculates prediction of abundances in n-member community by solving Ax=b where A is effective pairwise interactions and b set to -1
    abundance_n_member_community = slv(interactions_matrix, -1*np.ones((n,1)))

    return (interactions_matrix, abundance_n_member_community)