import numpy as np
import matplotlib.pyplot as plt

# Define the x and y ranges
x = np.linspace(-10, 10, 1000)
y = np.linspace(-10, 10, 1000)

# Create a meshgrid from x and y
X, Y = np.meshgrid(x, y)

# Define the five inequalities
ineq1 = X + 0 * Y <= 5.0
ineq2 = 1.3763784778814057 * X + Y <= 8.506491053308515
ineq3 = 0.32491882707222286 * X + Y <= 5.25730558133275
ineq4 = -0.32491882707222286 * X + Y <= 5.25730558133275
ineq5 = -1.3763784778814057 * X + Y <= 8.506491053308515

ineq6 = -1 * X + 0 * Y <= 5.0
ineq7 = -1.3763784778814057 * X - Y <= 8.506491053308515
ineq8 = -0.32491882707222286 * X - Y <= 5.25730558133275
ineq9 = 0.32491882707222286 * X - Y <= 5.25730558133275
ineq10 = 1.3763784778814057 * X - Y <= 8.506491053308515

# Find the boolean array of values that satisfy all five inequalities
intersection = (ineq1) & (ineq2) & (ineq3) & (ineq4) & (ineq5) & (ineq6) & (ineq7) & (ineq8) & (ineq9) & (ineq10)

# Define the plot
# plt.plot([], [], ' ', label='(1)*x+(0)*y<=5.0')
# plt.plot([], [], ' ', label='(0.32491882707222286)*x+(1)*y<=5.25730558133275')
# plt.plot([], [], ' ', label='(-1.3763784778814057)*x+(1)*y<=8.506491053308515')
# plt.plot([], [], ' ', label='(-1.3763784778814057)*x+(-1)*y<=8.506491053308515')
# plt.plot([], [], ' ', label='(0.32491882707222286)*x+(-1)*y<=5.25730558133275')
plt.contourf(X, Y, intersection, alpha=0.5)

# Add the legend and show the plot
plt.legend()
plt.show()
