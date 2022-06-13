# EPICS implement

Implementing the model described by Ansari et al in their article *An efficient and scalable top-down method for predicting structures of microbial communities*, NatCompSci 2021.

### How to
Number of species: Input n

Species abundance in monoculture: Input a vector array (MatLab) or list (Python)

Species abundances in leave-one-out cultures: 2D array where columns specify species and rows specify which leave-one-out setting it is

#### References
Model and n=8 data: Ansari et al 2021 *An efficient and scalable top-down method for predicting structures of microbial communities*, Nature Computational Science 1(9), pp. 619-628

n=5 data: Gould et al 2018 *Microbiome interactions shape host fitness*  Proceedings of the National Academy of Sciences (PNAS) 115(51), pp. E11951-E11960
