import numpy as np
from scipy.linalg import solve as slv

#Number of species
#number of species from Gould et al (as referenced by Ansari et al)
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
abundance_in_loo = [[0.0, 254186.34555361, 202495.347480754, 242561.433240061, 255231.590042967], [396525.520961999, 0.0, 205986.64657525, 272415.148100376, 241558.469147808], [123647.743095677, 21405.1659413566, 0.0, 59707.4297206303, 195981.399497278], [170548.611166451, 69566.7893094091, 104738.972834873, 0.0, 173192.864672013], [115120.312537355, 141809.224361488, 136160.664685335, 123146.5737988, 0.0]]
#below abundances from Ansari et al
#abundance_in_loo = [[0, 8317637.711, 1047128.548, 13.18256739, 281838293.1, 467735.1413, 20417.37945, 218776.1624], [91201.08394, 0, 48977.88194, 70794.57844, 27542.28703, 81283.05162, 89125.09381, 79432.82347], [16595869.07, 51286138.4, 0, 19952623.15, 1995262.315, 186208713.7, 407380.2778, 218776162.4], [831763.7711, 9549925.86, 35481338.92, 0, 3235936.569, 16595869.07, 50118.72336, 2630267.992], [1174.897555, 263026.7992, 95499.2586, 3162.27766, 0, 95499.2586, 22908.67653, 12882.49552], [36307.80548, 16595.86907, 46773.51413, 18620.87137, 288.4031503, 0, 123.0268771, 3019.95172], [4897.788194, 12022.64435, 23988.32919, 2630.267992, 4365.158322, 15848.93192, 0, 15848.93192], [11481.53621, 2290.867653, 331131.1215, 758.577575, 5623413.252, 630.9573445, 9549.92586, 0]]

#Creates matrix A for solving Ax=b, consists of species abundances in monoculture and l-o-o
abundance_matrix = np.zeros((n*n,n*n))
row = 0
for i in range(0,n):
    for j in range(0,n):
        if i==j:
            abundance_matrix[row][row] = abundance_in_monoc[i]
        else:
            abundance_matrix[row][n*j:n*j+n] = abundance_in_loo[i]
        row = row+1

#solves Ax=b creating a vector of effective pairwise interactions, b is a vector of -1 as growth rates are set to 1
interactions_vector = slv(abundance_matrix, -1*np.ones((n*n,1)))

#creates matrix of effective pairwise interactions to be used as A in Ax=b
interactions_matrix = np.reshape(interactions_vector,(n,n))

print("Interactions matrix:")
print(interactions_matrix)

#normalizes interaction matrix by self-interactions like done in EPICS paper
normalized_interactions_matrix = interactions_matrix
for i in range(0, n):
    for j in range(0, n):
        normalized_interactions_matrix[i][j] = normalized_interactions_matrix[i][j]/interactions_matrix[j][j]
print("\nNormalized interactions matrix:")
print(normalized_interactions_matrix)

#calculates prediction of abundances in n-member community by solving Ax=b where A is effective pairwise interactions and b set to -1
abundance_n_member_community = slv(interactions_matrix, -1*np.ones((n,1)))

print("\nPredicted abundances in an n-member community")
print(abundance_n_member_community)