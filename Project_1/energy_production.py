import numpy as np, matplotlib.pyplot as plt
import constants

"""
Temperature                 [T_9] = T_sol / 10^9 K  (unitless)
Reaction rate               [lambda_ik] = cm^3/s = 10^-6 m^3/s
Reaction rate per mass      [r_ik] = 1/kgs
Energy                      [Q_ik] = 10^-13 J = 10^-13 kgm^2/s^2
Energy production           [epsilon] = J/kgs = m^2/s^3

Upper limit of 7Be electron capture for T<10^6 is N_A lambda_e7 <= 1.57*10^-7/n_e (page 31 in lecture notes)
"""


print("Sanity check? (y/n): ", end="")
san_check = str(input())
#san_check = "n"
pp_check = 0
while san_check != "y" and san_check != "n":
    print("Not a valid input, try again! Sanity check? (y/n): ", end="")
    san_check = input()
print("")

def sanity_check(check):
    tol = 1e-2
    # Pre calculated values for sanity check
    if abs(check-4.04e2)/4.04e2 < tol:
        print("The check is valid for the initial pp reaction! Relative error: %.2f %%" % (100*abs(check-4.04e2)/4.04e2), "\n")
    elif abs(check-8.69e-9)/8.69e-9 < tol:
        print("The check is valid for the 33 reaction! Relative error: %.2f %%" % (100*abs(check-8.686e-9)/8.686e-9), "\n")
    elif abs(check-4.87e-5)/4.87e-5 < tol:
        print("The check is valid for the 34 reaction! Relative error: %.2f %%" % (100*abs(check-4.865e-5)/4.865e-5), "\n")
    elif abs(check-1.50e-6)/1.50e-6 < tol:
        print("The check is valid for the e7 reaction! Relative error: %.2f %%" % (100*abs(check-1.496e-6)/1.496e-6), "\n")
    elif abs(check-5.29e-4)/5.29e-4 < tol:
        print("The check is valid for the 17' reaction! Relative error: %.2f %%" % (100*abs(check-5.29e-4)/5.29e-4), "\n")
    elif abs(check-1.64e-6)/1.64e-6 < tol:
        print("The check is valid for the 17 reaction! Relative error: %.2f %%" % (100*abs(check-1.638e-6)/1.638e-6), "\n")
    elif abs(check-9.18e-8)/9.18e-8 < tol:
        print("The check is valid for the p14 reaction! Relative error: %.2f %%" % (100*abs(check-9.18e-8)/9.18e-8), "\n")
    else:
        print("The calculated value does not match any of the test results. The value is %.2e J m^-3 s^-1" % check, "\n")


def r_ik(X_i, X_k, A_i, A_k, lamb):
    # calculating the number densities for each of the colliding particles
    if X_i == "electron":
        n_i=C.n_etot
        n_k=C.rho_sol*X_k / (A_k*C.m_u)
    elif X_k == "electron":
        n_i=C.rho_sol*X_i / (A_i*C.m_u)
        n_k=C.n_etot
    else:
        n_i = C.rho_sol*X_i / (A_i*C.m_u)
        n_k = C.rho_sol*X_k / (A_k*C.m_u)
    # Kronecker-delta. If n_i=n_k, then the reaction rate should be divided by 2
    if n_i == n_k: dk = 1
    else: dk = 0
    return n_i*n_k*lamb / (C.rho_sol*(1+dk))


def pp():
    # initial reaction rates for each of the three pp-branches. Finds r*Q*rho
    r_pp = r_ik(C.X,C.X,1,1,C.l_pp)
    rQ_pp = r_pp*(C.Q_pp + C.Q_pd)*C.rho_sol

    # change the global value of pp_check
    global pp_check
    if san_check=="y" and pp_check==0:
        print("####### Initial pp reaction ########")
        print("pp sanity check:"); sanity_check(rQ_pp); pp_check=1
    return r_pp


def pp_b1():
    # reaction rate is gathered from initial pp reaction. Other
    # necessary reaction rates are included
    r_pp = pp()
    r_33 = r_ik(C.Y_3He,C.Y_3He,3,3,C.l_33)
    r_34 = r_ik(C.Y_3He,C.Y,3,4,C.l_34)

    # if the reaction rates of two 33 reaction + two 34 reaction is more
    # frequent than pp + pd reaction, 33 reaction is deacreased with the
    # ratio factor. Finds r*Q*rho
    if 2*r_33 + 2*r_34 > r_pp: r_33 = r_pp * r_33/(2*r_33+2*r_34)
    rQ_33 = r_33*C.Q_33*C.rho_sol
    if san_check=="y":
        print("####### pp I chain ########")
        print("33 sanity check:"); sanity_check(rQ_33)

    # generated energy is found for branch 1
    eps_b1 = rQ_33/C.rho_sol
    return eps_b1


def pp_b2():
    # reaction rate is gathered from initial pp reaction. Other
    # necessary reaction rates are included
    r_pp = pp()
    r_33 = r_ik(C.Y_3He,C.Y_3He,3,3,C.l_33)
    r_17 = r_ik(C.Z_7Be,C.X,7,1,C.l_17)

    # reaction rate test and calculating r*Q*rho
    r_34 = r_ik(C.Y_3He,C.Y,3,4,C.l_34)
    if 2*r_33 + 2*r_34 > r_pp: r_34 = r_pp * r_34/(2*r_33+2*r_34)
    rQ_34 = r_34*C.Q_34*C.rho_sol
    if san_check=="y":
        print("####### pp II chain ########")
        print("34 sanity check:"); sanity_check(rQ_34)

    # reaction e7 + 17 must be lower than 34 reaction as the last reaction
    # is necessary for both e7 and 17 reactions. Finds r*Q*rho
    r_e7 = r_ik(C.Z_7Be,"electron",7,0,C.l_e7)
    if r_e7+r_17 > r_34: r_e7 = r_e7*r_34/(r_e7+r_17)
    rQ_e7 = r_e7*C.Q_e7*C.rho_sol
    if san_check=="y": print("e7 sanity check:"); sanity_check(rQ_e7)

    # reaction rate test and calculating r*Q*rho
    r_17mrk = r_ik(C.Z_7Li,C.X,7,1,C.l_17mrk)
    if r_17mrk>r_e7: r_17mrk=r_e7
    rQ_17mrk = r_17mrk*C.Q_17mrk*C.rho_sol
    if san_check=="y": print("17' sanity check:"); sanity_check(rQ_17mrk)

    # generated energy is found for branch 2
    eps_b2 = (rQ_34 + rQ_e7 + rQ_17mrk)/C.rho_sol
    return eps_b2


def pp_b3():
    # reaction rate is gathered from initial pp reaction. Other
    # necessary reaction rates are included
    r_pp = pp()
    r_33 = r_ik(C.Y_3He,C.Y_3He,3,3,C.l_33)
    r_e7 = r_ik(C.Z_7Be,"electron",7,0,C.l_e7)

    # reaction rate test and calculating r*Q
    r_34 = r_ik(C.Y_3He,C.Y,3,4,C.l_34)
    if 2*r_33 + 2*r_34 > r_pp: r_34 = r_pp * r_34/(2*r_33+2*r_34)
    rQ_34 = r_34*C.Q_34*C.rho_sol
    if san_check=="y":
        print("####### pp III chain ########")
        print("34 sanity check:"); sanity_check(rQ_34)

    # reaction rate test and calculating r*Q*rho. Reaction e1 + 17 must
    # be lower than 34 reaction
    r_17 = r_ik(C.Z_7Be,C.X,7,1,C.l_17)
    if r_e7+r_17 > r_34: r_17 = r_17*r_34/(r_e7+r_17)

    # Now including rest of the Q's in this branch
    rQ_17 = r_17*(C.Q_17+C.Q_8+C.Q_8mrk)*C.rho_sol
    if san_check=="y": print("17 sanity check:"); sanity_check(rQ_17)

    # generated energy is found for branch 3
    eps_b3 = (rQ_34 + rQ_17)/C.rho_sol
    return eps_b3


def CNO():
    # finds reaction rate of p14 and calculating r*Q*rho
    r_p14 = r_ik(C.Z_14N,C.X,14,1,C.l_p14)
    rQ_p14 = r_p14*(C.Q_p12+C.Q_13+C.Q_p13+C.Q_p14+C.Q_15+C.Q_p15)*rho_sol
    if san_check=="y":
        print("####### CNO cycle chain ########")
        print("p14 sanity check:"); sanity_check(rQ_p14)

    # generated energy is found for CNO cycle
    eps_CNO = rQ_p14/rho_sol
    return eps_CNO



if __name__ == '__main__':
    # initialize Sun parameters
    T_sol = 1.57e7                      # K
    rho_sol = 1.62e5                    # kg/m^3
    C = constants.constants(T_sol,rho_sol)

    # energy generation for each branch and CNO cycle in the Sun
    b1 = pp_b1()
    b2 = pp_b2()
    b3 = pp_b3()
    CNO_ = CNO()
    tot = b1+b2+b3+CNO_
    print("Energy production in Sun pp I:        %.2e" % (b1/1.6022e-13), " MeV/kgs")
    print("Energy production in Sun pp II:       %.2e" % (b2/1.6022e-13), " MeV/kgs")
    print("Energy production in Sun pp III:      %.2e" % (b3/1.6022e-13), " MeV/kgs")
    print("Energy production in Sun CNO cycle:   %.2e" % (CNO_/1.6022e-13), " MeV/kgs \n")

    print("pp I relative energy production in Sun:        %.2f %%" % (100*b1/tot))
    print("pp II relative energy production in Sun:       %.2f %%" % (100*b2/tot))
    print("pp III relative energy production in Sun:      %.2f %%" % (100*b3/tot))
    print("CNO cycle relative energy production in Sun:   %.2f %%" % (100*CNO_/tot), "\n")

    # finds the total energy produced and amount of energy lost due to neutrinos
    Q_pp1 = (2*(C.Q_pp+C.Q_pd)+C.Q_33+2*C.Q_nu_pp)
    Q_pp2 = (C.Q_pp+C.Q_pd+C.Q_34+C.Q_e7+C.Q_17mrk+C.Q_nu_pp+C.Q_nu_e7)
    Q_pp3 = (C.Q_pp+C.Q_pd+C.Q_34+C.Q_17+C.Q_8+C.Q_8mrk+C.Q_nu_pp+C.Q_nu_8)
    Q_CNO = (C.Q_p12+C.Q_13+C.Q_p13+C.Q_p14+C.Q_15+C.Q_p15+C.Q_nu_13+C.Q_nu_15)
    print("Energy produced per chain reaction / lost energy due to neutrinos in pp I:        %.3f MeV / %.2f %%" % (Q_pp1/1.6022e-13, 100*2*C.Q_nu_pp/Q_pp1))
    print("Energy produced per chain reaction / lost energy due to neutrinos in pp II:       %.3f MeV / %.2f %%" % (Q_pp2/1.6022e-13, 100*(C.Q_nu_pp+C.Q_nu_e7)/Q_pp2))
    print("Energy produced per chain reaction / lost energy due to neutrinos in pp III:      %.3f MeV / %.2f %%" % (Q_pp3/1.6022e-13, 100*(C.Q_nu_pp+C.Q_nu_8)/Q_pp3))
    print("Energy produced per chain reaction / lost energy due to neutrinos in CNO cycle:   %.3f MeV / %.2f %%" % (Q_CNO/1.6022e-13, 100*(C.Q_nu_13+C.Q_nu_15)/Q_CNO))



    N = 3000
    T = np.logspace(4,9,N)
    eps_b1 = np.zeros(N)
    eps_b2 = np.zeros(N)
    eps_b3 = np.zeros(N)
    eps_CNO = np.zeros(N)
    eps_tot = 0


    # evolution of energy production at increasing temperature in core of the Sun
    san_check = "n"
    for i in range(N):
        C = constants.constants(T[i],rho_sol)
        eps_b1[i] = pp_b1()
        eps_b2[i] = pp_b2()
        eps_b3[i] = pp_b3()
        eps_CNO[i] = CNO()
        eps_tot = eps_b1[i]+eps_b2[i]+eps_b3[i]+eps_CNO[i]
        eps_b1[i] /= eps_tot
        eps_b2[i] /= eps_tot
        eps_b3[i] /= eps_tot
        eps_CNO[i] /= eps_tot


    plt.plot(T,eps_b1,label="pp I")
    plt.plot(T,eps_b2,label="pp II")
    plt.plot(T,eps_b3,label="pp III")
    plt.plot(T,eps_CNO,label="CNO cycle")
    plt.plot([T_sol,T_sol],[-0.05,1.05],"--",label="Sun")
    plt.title("Energy production in Sun as function of core temperature"); plt.legend(); plt.grid(); plt.xscale("log")
    plt.xlabel("Core temperature (K)"); plt.ylabel("Energy production $(\epsilon/\epsilon_{total})$")
    plt.show()
