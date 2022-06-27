from EPICS_function import epics
from GLV_steady_state_calculator import glv

(glv_interactions, glv_n_member_community, glv_monoc, glv_loo) = glv(5, 0.1, 0.2)
(epics_interactions, epics_normalized_interactions, epics_n_member_community) = epics(glv_monoc, glv_loo)

print("Interactions in GLV function")
print(glv_interactions)
print("\nInteractions calculated by EPICS")
print(epics_interactions)
print("\nInteractions used in GLV divided by interactions calculated in EPICS")
print(glv_interactions/epics_interactions)
print("\nInteractions used in GLV divided by normalized interactions calculated in EPICS")
print(glv_interactions/epics_normalized_interactions)
print("\nAbundances for n-member community from GLV")
print(glv_n_member_community)
print("\nAbundances for n-member community from epics")
print(epics_n_member_community)