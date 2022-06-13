# EPICS implement

Implementing the model described by Ansari et al in their article *An efficient and scalable top-down method for predicting structures of microbial communities*, NatCompSci 2021.

### How to use EPICS
Species abundance in monoculture: Input a vector array (MatLab) or list (Python) of species abundances in monoculture

Species abundances in leave-one-out cultures: 2D array where columns specify species and rows specify which leave-one-out setting it is.

### How to use generalized EPICS
Input is an excel-file with sheet 1 containing species abundance in monoculture in a row and sheet 2 containing species abundances in leave-one-out cultures, columns specify species and rows specify which leave-one-out setting it is.

### How to use GLV vs EPICS
Run code from index.py that calls a function for the generalized Lotka-Volterra model. The function simulates a community of n species with randomly assigned pairwise interactions between species. Simulated monoculture and leave-one-out culture abundances are used as an input for an EPICS function. Calculated pairwise interactions and abundances in n-member communities are printed out.

#### References
Model and n=8 data: Ansari et al 2021 *An efficient and scalable top-down method for predicting structures of microbial communities*, Nature Computational Science 1(9), pp. 619-628

n=5 data: Gould et al 2018 *Microbiome interactions shape host fitness*  Proceedings of the National Academy of Sciences (PNAS) 115(51), pp. E11951-E11960
