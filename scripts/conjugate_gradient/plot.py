
import matplotlib.pyplot as plt

# Data
gmg_v11 = [3.62e-02, 5.85e-03, 9.74e-04, 1.57e-04, 3.11e-05, 6.53e-06,
           1.43e-06, 3.25e-07, 7.52e-08, 1.78e-08, 4.26e-09, 1.04e-09,
           2.56e-10, 6.39e-11]
fmg1_v11 = [1.38e-06, 2.28e-08, 2.24e-09, 4.87e-10, 1.23e-10, 3.19e-11, 8.42e-12]
fmg2_v11 = [2.13e-08, 1.09e-09, 1.04e-10, 1.27e-11]

pcg_1xf = [3.62219e-02, 9.3984e-07, 8.64942e-09, 5.94761e-10, 3.95338e-11]
fmg1_pcg = [1.37916e-06, 2.03839e-08, 1.09761e-09, 1.02575e-10, 8.13779e-12]
fmg2_pcg = [2.13099e-08, 8.93099e-10, 2.99579e-11]

# Iteration indices
iters_gmg_v11 = list(range(len(gmg_v11)))
iters_pcg_1xf = list(range(len(pcg_1xf)))
iters_fmg1_v11 = list(range(len(fmg1_v11)))
iters_fmg1_pcg = list(range(len(fmg1_pcg)))
iters_fmg2_v11 = list(range(len(fmg2_v11)))
iters_fmg2_pcg = list(range(len(fmg2_pcg)))

# Plot
plt.figure(figsize=(8,6))

# GMGPolar (solid lines)
plt.semilogy(iters_gmg_v11, gmg_v11, '-o', color='tab:blue', label="GMGPolar. No FMG, Cycle: V(1,1)")
plt.semilogy(iters_pcg_1xf, pcg_1xf, '--o', color='tab:blue', label="MG-PCG. No FMG, Cycle: 1x F-FMG + No MG-Iter.")
plt.semilogy(iters_fmg1_v11, fmg1_v11, '-s', color='tab:green', label="GMGPolar. 1x F-FMG, Cycle: V(1,1)")
plt.semilogy(iters_fmg1_pcg, fmg1_pcg, '--s', color='tab:green', label="MG-PCG. 1x F-FMG, Cycle: 1x F-FMG + No MG-Iter.")
plt.semilogy(iters_fmg2_v11, fmg2_v11, '-^', color='tab:orange', label="GMGPolar. 2x F-FMG, Cycle: V(1,1)")
plt.semilogy(iters_fmg2_pcg, fmg2_pcg, '--^', color='tab:orange', label="MG-PCG. 2x F-FMG, Cycle: 1x F-FMG + No MG-Iter.")

# Formatting
plt.xlabel("Iteration")
plt.ylabel(r"$\|r_k\|$")
plt.title("Residual Convergence of GMGPolar vs GMGPolar-PCG\n" "Config: Czarny, CartesianR6, No Extrapolation (1537 × 2048)")
plt.legend()
plt.grid(True, which="both", ls="--", alpha=0.6)

# Save to file
plt.tight_layout()
plt.savefig("gmgpolar_convergence.png", dpi=300)
plt.show()
