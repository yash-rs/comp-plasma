import csv
import numpy as np
import matplotlib.pyplot as plt

h, l1, l2, linf = [], [], [], []

with open("data/poisson_convergence.csv") as f:
    for row in csv.DictReader(f):
        h.append(float(row["h"]))
        l1.append(float(row["l1"]))
        l2.append(float(row["l2"]))
        linf.append(float(row["linf"]))

h = np.array(h)

for name, err in [("l1", l1), ("l2", l2), ("linf", linf)]:
    err = np.array(err)
    slope, intercept = np.polyfit(np.log(h), np.log(err), 1)
    print(f"{name}: slope = {slope:.3f}")
    plt.loglog(h, err, "o-", label=f"{name} (slope={slope:.2f})")

plt.xlabel("h")
plt.ylabel("error")
plt.title("Convergence of FD Poisson solver")
plt.legend()
plt.savefig("data/poisson_convergence.png")
