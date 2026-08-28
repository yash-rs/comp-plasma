import csv
import matplotlib.pyplot as plt

x, phi_numeric, phi_exact = [], [], []

with open("data/poisson_solution.csv") as f:
    for row in csv.DictReader(f):
        x.append(float(row["x"]))
        phi_numeric.append(float(row["phi_numeric"]))
        phi_exact.append(float(row["phi_exact"]))

plt.plot(x, phi_exact, label="analytical")
plt.plot(x, phi_numeric, "--", label="numerical")
plt.xlabel("x")
plt.ylabel("phi(x)")
plt.title("1D Poisson solution")
plt.legend()
plt.savefig("data/poisson_solution.png")
