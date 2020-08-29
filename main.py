'''
This is working!! now just incorporate the EC fit correctly
'''
from AnharmonicOscillator import *

E0 = AnharmonicOscillator(nEigState = 1) # ground state, lowest; nEigvectors is howmany eigvects to use for EC 
E1 = AnharmonicOscillator(nEigState = 2) # 2nd lowest stae
E3 = AnharmonicOscillator(nEigState = 3) # 3rd lowest state
E4 = AnharmonicOscillator(nEigState = 4)
#the calculated points
figure("3 lowest eigen states vs. λ")
scatter(E0.λ_list, E0.eigValues, label="$E_0$ (calculated)", marker="s") # ploting the ground state
scatter(E1.λ_list, E1.eigValues, label="$E_1$ (calculated)", marker="s") # ploting the 2nd state
scatter(E3.λ_list, E3.eigValues, label="$E_3$ (calculated)", marker="s")
scatter(E4.λ_list, E4.eigValues, label="$E_4$ (calculated)", marker="s")
title("Eigen States vs. $\\lambda$")
xlabel("$\\lambda$")
ylabel("Eigen State")
legend()
show()

# showing the tabels
E0.prettyPrint()
E1.prettyPrint()
E3.prettyPrint()
E4.prettyPrint()