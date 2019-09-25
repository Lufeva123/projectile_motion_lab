import pandas as pd
import matplotlib.pyplot as plt
import math

################################################################################################################################################################
# Define functions to calculate x
g = 9.807
DECIMAL_PLACES = 4

def calculate_delta_h (h1, h2):
    return round(h2 - h1, DECIMAL_PLACES)


def calculate_delta_h_prime(h1_prime, h2_prime):
    return round(h2_prime - h1_prime, DECIMAL_PLACES)


def calculate_v_0 (delta_h, delta_h_prime):
    return round(math.sqrt(10/7 * g * (delta_h_prime - delta_h)), DECIMAL_PLACES)


def calculate_v_0_x(v_0, d, l):
    return round((v_0 * d)/l, DECIMAL_PLACES)


def calculate_v_0_y(v_0, h2, h3, l):
    return round((v_0 * (h2 - h3))/l, DECIMAL_PLACES)


def calculate_t(v_0_y, h2):
    delta = (v_0_y ** 2) + (2 * g * h2)
    return round((v_0_y + math.sqrt(delta))/g, DECIMAL_PLACES)

def calculate_range(v_0_x, t):
    return round(v_0_x * t, DECIMAL_PLACES)

################################################################################################################################################################
# Define functions to calculate uncertainty in x
def calculate_x(h1_v, h1_prime_v, h2_v, h2_prime_v, h3_v, d_v, l_v):
    # Unpacking values and uncertainties
    h1, u_h1 = h1_v
    h1_prime, u_h1_prime = h1_prime_v
    h2, u_h2 = h2_v
    h2_prime, u_h2_prime = h2_prime_v
    h3, u_h3 = h3_v
    d, u_d = d_v
    l, u_l = l_v

    # Calculating delta_h and delta_h_prime with uncertainties (u)
    delta_h = calculate_delta_h(h1, h2)
    u_delta_h = math.sqrt((u_h1**2) + (u_h2**2))
    delta_h_prime = calculate_delta_h_prime(h1_prime, h2_prime)
    u_delta_h_prime = math.sqrt((u_h1_prime**2) + (u_h2_prime**2))

    # Calculating v_0 with uncertainty (u)
    v_0 = calculate_v_0(delta_h, delta_h_prime)
    deltas_subtraction = delta_h_prime - delta_h
    u_deltas_subtraction = math.sqrt((u_delta_h**2) + (u_delta_h_prime**2))
    u_v_0 = v_0 * (1/2) * (u_deltas_subtraction/deltas_subtraction)

    # Calculating v_0_x and uncertainty (u)
    v_0_x = calculate_v_0_x(v_0, d, l)
    u_v_0_x = v_0_x * ((u_v_0/v_0) + (u_d/d) + (u_l/l))

    #Caclulating v_0_y and uncertainty (u)
    v_0_y = calculate_v_0_y(v_0, h2, h3, l)
    u_h2_minus_h3 = math.sqrt((u_h2**2) + (u_h3**2))
    u_v_0_y = v_0_y * ((u_v_0/v_0) + (u_h2_minus_h3/abs(h2 - h3)) + (u_l/l))

    # Caclulating t and uncertainty (u)
    t = calculate_t(v_0_y, h2)
    delta = (v_0_y ** 2) + (2 * g * h2)
    u_v_0_y_squared = (v_0_y**2) * 2 * (u_v_0_y/v_0_y)
    u_delta = math.sqrt((u_v_0_y_squared**2) + (2 * g * u_h2))
    u_sqrt_delta = u_delta * (1/2) * (u_delta/delta)
    u_t = math.sqrt((u_v_0_y**2) + (u_sqrt_delta**2))/g

    # Calculating range and uncertainty (u)
    rng = calculate_range(v_0_x, t)
    u_rng = rng * ((u_v_0_x/v_0_x) + (u_t/t))

    return (rng, u_rng)


################################################################################################################################################################
# Read the lab data in as data frames
heavy_trials_1_df = pd.read_csv("data/lab_1_ heavy _trials_1.csv")
heavy_trials_2_df = pd.read_csv("data/lab_1_ heavy _trials_2.csv")
plastic_trials_1_df = pd.read_csv("data/lab_1_ plastic _trials_1.csv")
plastic_trials_2_df = pd.read_csv("data/lab_1_ plastic _trials_2.csv")

################################################################################################################################################################
# Calculate uncertainties of expected horizontal landing position x
h1_u = 0.2 #cm
h1_prime_u = 0.2 #cm
h2_u = 0.2 #cm
h2_prime_u = 0.2 #cm
h3_u = 0.2 #cm
d_u = 0.1 #cm
l_u = 0.1 #cm

# Heavy ball trials 1
print("Heavy ball trial 1")
h1, h1_prime, h2, h2_prime, h3, d, l = heavy_trials_1_df[["h1 (cm)", "h'1 (cm)", "h2 (cm)", "h'2 (cm)", "h3 (cm)", "D (cm)", "L (cm)"]].iloc[0]
heavy_1_x_pred, heavy_1_x_pred_u = calculate_x((h1, h1_u), (h1_prime, h1_prime_u), (h2, h2_u), (h2_prime, h2_prime_u), (h3, h3_u), (d, d_u), (l, l_u))
heavy_1_x_pred, heavy_1_x_pred_u = round(heavy_1_x_pred, 1), round(heavy_1_x_pred_u, 1)
print("     Predicted x:", heavy_1_x_pred, "/ Uncertainty in x:", heavy_1_x_pred_u)

# Heavy ball trials 2
print("\nHeavy ball trial 2")
h1, h1_prime, h2, h2_prime, h3, d, l = heavy_trials_2_df[["h1 (cm)", "h'1 (cm)", "h2 (cm)", "h'2 (cm)", "h3 (cm)", "D (cm)", "L (cm)"]].iloc[0]
heavy_2_x_pred, heavy_2_x_pred_u = calculate_x((h1, h1_u), (h1_prime, h1_prime_u), (h2, h2_u), (h2_prime, h2_prime_u), (h3, h3_u), (d, d_u), (l, l_u))
heavy_2_x_pred, heavy_2_x_pred_u = round(heavy_2_x_pred, 1), round(heavy_2_x_pred_u, 1)
print("     Predicted x:", heavy_2_x_pred, "/ Uncertainty in x:", heavy_2_x_pred_u)

# Plastic ball trials 1
print("\nPlastic ball trials 1")
h1, h1_prime, h2, h2_prime, h3, d, l = plastic_trials_1_df[["h1 (cm)", "h'1 (cm)", "h2 (cm)", "h'2 (cm)", "h3 (cm)", "D (cm)", "L (cm)"]].iloc[0]
plastic_1_x_pred, plastic_1_x_pred_u = calculate_x((h1, h1_u), (h1_prime, h1_prime_u), (h2, h2_u), (h2_prime, h2_prime_u), (h3, h3_u), (d, d_u), (l, l_u))
plastic_1_x_pred, plastic_1_x_pred_u = round(plastic_1_x_pred, 1), round(plastic_1_x_pred_u, 1)
print("     Predicted x:", plastic_1_x_pred, "/ Uncertainty in x:", plastic_1_x_pred_u)

# Plastic ball trials 2
print("\nPlastic ball trials 2")
h1, h1_prime, h2, h2_prime, h3, d, l = plastic_trials_2_df[["h1 (cm)", "h'1 (cm)", "h2 (cm)", "h'2 (cm)", "h3 (cm)", "D (cm)", "L (cm)"]].iloc[0]
plastic_2_x_pred, plastic_2_x_pred_u = calculate_x((h1, h1_u), (h1_prime, h1_prime_u), (h2, h2_u), (h2_prime, h2_prime_u), (h3, h3_u), (d, d_u), (l, l_u))
plastic_2_x_pred, plastic_2_x_pred_u = round(plastic_2_x_pred, 1), round(plastic_2_x_pred_u, 1)
print("     Predicted x:", plastic_2_x_pred, "/ Uncertainty in x:", plastic_2_x_pred_u)

################################################################################################################################################################
# Calculate mean positions and standard error of mean for x and z in each trial

# Heavy ball trials 1
print("\n\nHeavy ball trial 1")
heavy_1_x_values = heavy_trials_1_df["x (actual) (cm)"] + heavy_1_x_pred
heavy_1_x_mean = heavy_1_x_values.mean()
heavy_1_x_std = heavy_1_x_values.std()
heavy_1_x_sem = heavy_1_x_values.sem()
heavy_1_x_mean, heavy_1_x_sem = round(heavy_1_x_mean, 1), round(heavy_1_x_sem, 1)
print("     Mean x:", heavy_1_x_mean, "/ Standard error of the mean:", heavy_1_x_sem)

heavy_1_z_values = heavy_trials_1_df["z (cm)"]
heavy_1_z_mean = heavy_1_z_values.mean()
heavy_1_z_std = heavy_1_z_values.std()
heavy_1_z_sem = heavy_1_z_values.sem()
heavy_1_z_mean, heavy_1_z_sem = round(heavy_1_z_mean, 1), round(heavy_1_z_sem, 1)
print("     Mean z:", heavy_1_z_mean, "/ Standard error of the mean:", heavy_1_z_sem)
print(heavy_1_x_std, heavy_1_z_std)

# Heavy ball trials 2
print("\nHeavy ball trial 2")
heavy_2_x_values = heavy_trials_2_df["x (actual) (cm)"] + heavy_2_x_pred
heavy_2_x_mean = heavy_2_x_values.mean()
heavy_2_x_std = heavy_2_x_values.std()
heavy_2_x_sem = heavy_2_x_values.sem()
heavy_2_x_mean, heavy_2_x_sem = round(heavy_2_x_mean, 1), round(heavy_2_x_sem, 1)
print("     Mean x:", heavy_2_x_mean, "Standard error of the mean:", heavy_2_x_sem)

heavy_2_z_values = heavy_trials_2_df["z (cm)"]
heavy_2_z_mean = heavy_2_z_values.mean()
heavy_2_z_std = heavy_2_z_values.std()
heavy_2_z_sem = heavy_2_z_values.sem()
heavy_2_z_mean, heavy_2_z_sem = round(heavy_2_z_mean, 1), round(heavy_2_z_sem, 1)
print("     Mean z:", heavy_2_z_mean, "/ Standard error of the mean:", heavy_2_z_sem)
print(heavy_2_x_std, heavy_2_z_std)

# Plastic ball trials 1
print("\nPlastic ball trial 1")
plastic_1_x_values = plastic_trials_1_df["x (actual) (cm)"] + plastic_1_x_pred
plastic_1_x_mean = plastic_1_x_values.mean()
plastic_1_x_std = plastic_1_x_values.std()
plastic_1_x_sem = plastic_1_x_values.sem()
plastic_1_x_mean, plastic_1_x_sem = round(plastic_1_x_mean, 1), round(plastic_1_x_sem, 1)
print("     Mean x:", plastic_1_x_mean, "/ Standard error of the mean:", plastic_1_x_sem)

plastic_1_z_values = plastic_trials_1_df["z (cm)"]
plastic_1_z_mean = plastic_1_z_values.mean()
plastic_1_z_std = plastic_1_z_values.std()
plastic_1_z_sem = plastic_1_z_values.sem()
plastic_1_z_mean, plastic_1_z_sem = round(plastic_1_z_mean, 1), round(plastic_1_z_sem, 1)
print("     Mean z:", plastic_1_z_mean, "/ Standard error of the mean:", plastic_1_z_sem)
print(plastic_1_x_std, plastic_1_z_std)
# Plastic ball trials 2
print("\nPlastic ball trial 2")
plastic_2_x_values = plastic_trials_2_df["x (actual) (cm)"] + plastic_2_x_pred
plastic_2_x_mean = plastic_2_x_values.mean()
plastic_2_x_std = plastic_2_x_values.std()
plastic_2_x_sem = plastic_2_x_values.sem()
plastic_2_x_mean, plastic_2_x_sem = round(plastic_2_x_mean, 1), round(plastic_2_x_sem, 1)
print("     Mean x:", plastic_2_x_mean, "/ Standard error of the mean:", plastic_2_x_sem)

plastic_2_z_values = plastic_trials_2_df["z (cm)"]
plastic_2_z_mean = plastic_2_z_values.mean()
plastic_2_z_std = plastic_2_z_values.std()
plastic_2_z_sem = plastic_2_z_values.sem()
plastic_2_z_mean, plastic_2_z_sem = round(plastic_2_z_mean, 1), round(plastic_2_z_sem, 1)
print("     Mean z:", plastic_2_z_mean, "/ Standard error of the mean:", plastic_2_z_sem)
print(plastic_2_x_std/math.sqrt(20), plastic_2_z_std)
################################################################################################################################################################
# Create histogram of the difference between predicted x and actual x

# Heavy ball 1 trials histogram
# x histogram
heavy_1_x_figure = plt.figure()
heavy_1_x_displacement = heavy_trials_1_df["x (actual) (cm)"]
plt.hist(heavy_1_x_displacement)
plt.title("Distribution of x displacements for heavy ball trial 1")
plt.xlabel("Actual range displacement from predicted range (cm)")
plt.ylabel("Frequency")
plt.savefig("graphs/heavy_1_x_displacement_histogram.png")
# z histogram
heavy_1_z_figure = plt.figure()
heavy_1_z_displacement = heavy_trials_1_df["z (cm)"]
plt.hist(heavy_1_z_displacement)
plt.title("Distribution of z displacements for heavy ball trial 1")
plt.xlabel("Actual range displacement from predicted range (cm)")
plt.ylabel("Frequency")
plt.savefig("graphs/heavy_1_z_displacement_histogram.png")

# Heavy ball 2 trials histogram
# x histgogram
heavy_2_x_figure = plt.figure()
heavy_2_x_displacement = heavy_trials_2_df["x (actual) (cm)"]
plt.hist(heavy_2_x_displacement)
plt.title("Distribution of x displacements for heavy ball trial 2")
plt.xlabel("Actual range displacement from predicted range (cm)")
plt.ylabel("Frequency")
plt.savefig("graphs/heavy_2_x_displacement_histogram.png")
# z histogram
heavy_2_z_figure = plt.figure()
heavy_2_z_displacement = heavy_trials_2_df["z (cm)"]
plt.hist(heavy_2_z_displacement)
plt.title("Distribution of z displacements for heavy ball trial 2")
plt.xlabel("Actual range displacement from predicted range (cm)")
plt.ylabel("Frequency")
plt.savefig("graphs/heavy_2_z_displacement_histogram.png")

# Plastic ball 1 trials histogram
# x histogram
plastic_1_x_figure = plt.figure()
plastic_1_x_displacement = plastic_trials_1_df["x (actual) (cm)"]
plt.hist(plastic_1_x_displacement)
plt.title("Distribution of x displacements for plastic ball trial 1")
plt.xlabel("Actual range displacement from predicted range (cm)")
plt.ylabel("Frequency")
plt.savefig("graphs/plastic_1_x_displacement_histogram.png")
# z histogram
plastic_1_z_figure = plt.figure()
plastic_1_z_displacement = plastic_trials_1_df["z (cm)"]
plt.hist(plastic_1_z_displacement)
plt.title("Distribution of z displacements for plastic ball trial 1")
plt.xlabel("Actual range displacement from predicted range (cm)")
plt.ylabel("Frequency")
plt.savefig("graphs/plastic_1_z_displacement_histogram.png")

# Plastic ball 2 trials histogram
# x histogram
plastic_2_x_figure = plt.figure()
plastic_2_x_displacement = plastic_trials_2_df["x (actual) (cm)"]
plt.hist(plastic_2_x_displacement)
plt.title("Distribution of x displacements for plastic ball trial 2")
plt.xlabel("Actual range displacement from predicted range (cm)")
plt.ylabel("Frequency")
plt.savefig("graphs/plastic_2_x_displacement_histogram.png")
# z histogram
plastic_2_z_figure = plt.figure()
plastic_2_z_displacement = plastic_trials_2_df["z (cm)"]
plt.hist(plastic_2_z_displacement)
plt.title("Distribution of z displacements for plastic ball trial 2")
plt.xlabel("Actual range displacement from predicted range (cm)")
plt.ylabel("Frequency")
plt.savefig("graphs/plastic_2_z_displacement_histogram.png")
