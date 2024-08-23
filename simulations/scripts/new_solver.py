import numpy as np
from scipy.special import gamma
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import minimize

def solve_equilibrium(T_c, T_k, t_ck, theta, rho):
    # Set the initial values and flatten for solver
    w_ck = np.ones((10,5)).flatten()

    def objective(w_ck):
        # Reshape the input to the correct shape
        w_ck = w_ck.reshape((10, 5))
        # Define Q and Y
        Q = np.outer(T_c, T_k) * gamma((theta - 1) / theta) * t_ck * w_ck ** (theta - 1) * np.sum(t_ck * w_ck ** theta, axis=0) ** (-rho)
        Q = np.sum(Q)
        Y = np.sum(gamma((theta - 1) / theta) * t_ck * w_ck ** theta * np.sum(t_ck * w_ck ** theta, axis=0) ** (-rho), axis=1) / np.sum(w_ck / np.outer(T_c, T_k), axis=1)
        Y = np.sum(Y)
        # Return the objective function which is the difference between Q and Y
        return (Q - Y) ** 2

    # Minimize the objective function
    w_ck = minimize(objective, w_ck, method='BFGS').x.reshape((10, 5))
    w_ck = np.abs(w_ck)

    # Calculate the rest of the equilibrium
    lambda_k = np.sum(t_ck * w_ck, axis=0)
    pi_ck_left = t_ck * (w_ck ** theta)
    pi_ck_right = (lambda_k ** (- rho)) / np.sum(lambda_k ** (1 - rho))
    pi_ck = pi_ck_left * pi_ck_right

    pi_c_left = (t_ck * (w_ck ** theta)) / np.sum(t_ck * (w_ck ** theta), axis=0)
    pi_c_right = (lambda_k ** (1 - rho)) / np.sum(lambda_k ** (1 - rho))
    pi_c = np.sum(pi_c_left * pi_c_right, axis=1)

    return pi_c, pi_ck

np.random.seed(42)

T_c = np.random.normal(loc=10, scale=3, size=10)
T_k = np.random.normal(loc=10, scale=3, size=5)
t_ck = np.random.normal(loc=10, scale=3, size=(10,5))

similarity_factor = 0.8  # Factor to control similarity (1.0 means identical)

for i in [1, 3, 8]:
    t_ck[i] = similarity_factor * t_ck[6] + (1 - similarity_factor) * np.random.normal(loc=10, scale=3, size=(5,))

theta = 3
rho = 0.4
alpha = 7
eta = 1.65

# Calculate the equilibrium using the new solver
pi_c, pi_ck = solve_equilibrium(T_c, T_k, t_ck, theta, rho)

T_c[6] = T_c[6] * 0.85
pi_shock_c, pi_shock_ck = solve_equilibrium(T_c, T_k, t_ck, theta, rho)

reference_row = t_ck[6]
correlation = []
# Calculate correlation excluding row 6
for i, row in enumerate(t_ck):
    if i != 6:
        correlation.append(np.abs(np.corrcoef(reference_row, row)[0, 1]))

# Calculate change
change = pi_shock_c / pi_c

# Ensure change has the same length as correlation (excluding row 6)
change = np.delete(change, 6)

sizes = np.delete(pi_c, 6)

sns.set(style="whitegrid")
plt.figure(figsize=(10, 6))

# Assuming correlation and change are already defined
# Plot correlation vs. change using seaborn
sns.scatterplot(x=correlation, y=change, size=sizes, sizes=(20, 200))
plt.xlabel('Correlation')
plt.ylabel('Change')
plt.title('Correlation vs. Change')
plt.legend().set_visible(False)
# plt.savefig('graphs/city_shock.png')
plt.show()

print(f'pi_c: {pi_c}')
print(f'pi_shock_c: {pi_shock_c}')