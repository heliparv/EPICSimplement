import numpy as np
from scipy.linalg import solve as slv
from scipy.stats import bernoulli

#Function finds a steady state for a generalized Lotka-Volterra model with only pairwise interactions.
#Number of species is an input value.
#Sets random pairwise coefficients from a normal distribution, standard deviation is an input value of the function
#Introduces sparcity by setting some interactions to zero, input value sparcity defines how much sparcity is introduced
#returns interactions matrix, abundance in n-member community, monoculture abundances and leave-one-out abundances

def glv(species, std, sparcity):
    #Number of species
    n = species

    #intrinsic growth rate for species (same for all species)
    ri = 1/n

    #create interactions matrix by drwing from normal distribution
    interactions_matrix = np.random.normal(loc=0.0, scale=std, size=(n,n))

    #introduce sparcity by setting some interactions to zero with boolean indexing
    draw = bernoulli(sparcity)
    sparc_matr = np.array(np.reshape(draw.rvs(n*n), (n, n)), dtype=bool)
    interactions_matrix[sparc_matr] = 0

    #self-interaction coefficients are -1
    for i in range(0,n):
        interactions_matrix[i][i] = -1

    #calculate steady state abundances in n-member community
    abundance_n_member_community = slv(interactions_matrix, -ri*np.ones((n,1)))

    #calculate steady state abundances in leave-one-out cultures
    loo_abundances = np.zeros((n,n))
    for i in range(0, n):
        temp_interactions = np.delete(interactions_matrix, i, 0)
        temp_interactions = np.delete(temp_interactions, i, 1)
        abund = slv(temp_interactions, -ri*np.ones((n-1,1)))
        abund = np.insert(abund, i, 0)
        loo_abundances[i] = abund

    #monoculture abundances are all 1/n due to assumptions ri=1/n and self-interaction = -1
    monoc_abundances = (1/n)*np.ones(n)

    return (interactions_matrix, abundance_n_member_community, monoc_abundances, loo_abundances)