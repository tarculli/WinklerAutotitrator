import winkler_functions as wf
import numpy as np  # type: ignore
from scipy.optimize import curve_fit  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import pandas as pd  # type: ignore

filepath = input(f"What is the filepath for the data you want to plot?\n")

# Read the CSV file while skipping comment lines
try:
    data = pd.read_csv(filepath, comment='#')  # Skip metadata lines
except Exception as e:
    print(f"Error reading CSV: {e}")
    data = None

# Proceed only if data was read successfully
if data is not None:
    # Extract data
    volume_data = data["Volume (mL)"]
    potential_data = data["Potential (mV)"]

    # Step 1: Calculate absolute value of the derivative of the potential
    dy_dx = np.abs((np.nan_to_num(np.gradient(potential_data, volume_data), nan=0)))

    # Step 2: Define the Cauchy (Lorentzian) function
    def cauchy(x, A, x0, gamma):
        return A / (1 + ((x - x0) / gamma) ** 2)

    # Step 3: Fit the Cauchy function to the absolute value of the derivative
    initial_guess = [max(dy_dx), volume_data[np.argmax(dy_dx)], 0.001]

    try:
        # Enforce A >= max of derivative, gamma > 0, x0 within volume range
        lower_bounds = [max(dy_dx), min(volume_data), 1e-6]
        upper_bounds = [np.inf, max(volume_data), np.inf]

        params, covariance = curve_fit(
            cauchy,
            volume_data,
            dy_dx,
            p0=initial_guess,
            bounds=(lower_bounds, upper_bounds)
        )

        A_fitted, x0_fitted, gamma_fitted = params
        perr = np.sqrt(np.diag(covariance))  # Parameter standard deviations
        A_err, x0_err, gamma_err = perr
        fitted_curve = cauchy(np.linspace(0, max(volume_data), 10000), *params)

    except Exception as e:
        print(f"Curve fitting failed: {e}")
        params = None


    # Create subplots
    fig, ax = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

    # Top Panel: Raw Titration Data
    ax[0].scatter(volume_data, potential_data, color='blue', label="Raw Data", s=12)
    ax[0].set_ylabel("Potential (mV)")
    ax[0].set_title(f"{filepath[5:-4]}")

    # Add a line for the endpoint location
    if params is not None:
        ax[0].axvline(x0_fitted, color='green', linestyle=':', label=f"Endpoint (x₀) = {x0_fitted:.3f} mL")

    ax[0].legend()
    ax[0].grid(True)

    # Bottom Panel: Derivative and Cauchy Fit
    ax[1].scatter(volume_data, dy_dx, label="First Derivative", color='blue', s=12)

    if params is not None:
        ax[1].plot(np.linspace(0, max(volume_data), 10000), fitted_curve, label="Lorentzian Fit", color='red', linestyle='--')
        ax[1].axvline(x0_fitted, color='green', linestyle=':', label=f"Endpoint (x₀) = {x0_fitted:.3f} ± {2*x0_err:.3f} mL")

        # Display fit results in a textbox
        fit_text = (
            f"Fit Parameters:\n"
            f"A = {A_fitted:.1f} ± {A_err:.1f}\n"
            f"x₀ = {x0_fitted:.3f} ± {x0_err:.3f} mL\n"
            f"γ = {gamma_fitted:.3f} ± {gamma_err:.3f} mL"
        )
        ax[1].text(
            0.25, 0.95, fit_text,
            transform=ax[1].transAxes,
            fontsize=10,
            verticalalignment='top',
            horizontalalignment='right',
            bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.8)
        )

    ax[1].set_xlabel("Titrant Volume (mL)")
    ax[1].set_ylabel(r"Absolute Value of Derivative")
    ax[1].set_title("Lorentzian Fit to Derivative of Titration Data")
    ax[1].legend()
    ax[1].grid(True)

    plt.tight_layout()
    plt.show()
