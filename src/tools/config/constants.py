############## CONSTANTS ##############

G = 6.674e-11          # m^3 / (kg s^2)
C = 3.0e8              # m/s
C_SQ = C ** 2
G_OVER_C2 = G / C_SQ

KPC_TO_M = 3.086e19
MPC_TO_M = 3.086e22
M_SUN = 1.989e30

GYR_TO_S = 1e9 * 365.25 * 24 * 3600
HUBBLE_H = 0.7         # your assumed h
R0_KPC = 0.5           # galaxy softening radius
R0_M = R0_KPC * KPC_TO_M

KM_S_TO_M_S = 1000.0