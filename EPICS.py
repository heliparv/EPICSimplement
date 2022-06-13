import numpy as np
from scipy.linalg import solve as slv
from copy import deepcopy

#Number of species

#below number of species from Gould et al (as referenced by Ansari et al)
n = 5

#below number of species in study by Ansari et al
#n=8


#abundance in monoculture

#below abundances from Gould et al (as referenced by Ansari et al)
abundance_in_monoc = [177827.941003893, 422668.614265604, 143301.257023696, 237137.370566166, 143301.257023696]

#below abundances from Ansari et al
#abundance_in_monoc = [2884030259.24548, 41686.9890738402, 1258925.78380719, 436515.384984744, 8709663.37151069, 199526.324505624, 346737.031168192, 2137962.73530952]


#abundances in leave-one-out cultures

#below abundances from Gould et al (as referenced by Ansari et al)
abundance_in_loo = np.array([[0.0, 396525.520961999, 123647.743095677, 170548.611166451, 115120.312537355], [254186.34555361, 0.0, 21405.1659413566, 69566.7893094091, 141809.224361488], [202495.347480754, 205986.64657525, 0.0, 104738.972834873, 136160.664685335], [242561.433240061, 272415.148100376, 59707.4297206303, 0.0,  123146.5737988], [255231.590042967,  241558.469147808,  195981.399497278,  173192.864672013, 0.0]])

#below abundances from Ansari et al
#abundance_in_loo = [[0, 91201.08394, 16595869.07, 831763.7711, 1174.897555, 36307.80548, 4897.788194, 11481.53621], [8317637.711, 0, 51286138.4, 9549925.86, 263026.7992, 16595.86907, 12022.64435, 2290.867653], [1047128.548, 48977.88194, 0, 35481338.92, 95499.2586, 46773.51413, 23988.32919, 331131.1215], [13.18256739, 70794.57844, 19952623.15, 0, 3162.27766, 18620.87137, 2630.267992, 758.577575], [281838293.1, 27542.28703, 1995262.315, 3235936.569, 0, 288.4031503, 4365.158322, 5623413.252], [467735.1413, 81283.05162, 186208713.7, 16595869.07, 95499.2586, 0, 15848.93192, 630.9573445], [20417.37945, 89125.09381, 407380.2778, 50118.72336, 22908.67653, 123.0268771, 0, 9549.92586], [218776.1624, 79432.82347, 218776162.4, 2630267.992, 12882.49552, 3019.95172, 15848.93192, 0]]

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

print("Interactions matrix:")
print(interactions_matrix)

#normalizes interaction matrix by self-interactions like done in EPICS paper
normalized_interactions_matrix = deepcopy(interactions_matrix)
for i in range(0, n):
    for j in range(0, n):
        normalized_interactions_matrix[i][j] = normalized_interactions_matrix[i][j]/interactions_matrix[j][j]
print("\nNormalized interactions matrix:")
print(normalized_interactions_matrix)

#calculates prediction of abundances in n-member community by solving Ax=b where A is effective pairwise interactions and b set to -1
abundance_n_member_community = slv(interactions_matrix, -1*np.ones((n,1)))

print("\nPredicted abundances in an n-member community")
print(abundance_n_member_community)