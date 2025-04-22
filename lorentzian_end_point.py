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

if data is not None:
    # Average potentials for repeated volume values
    grouped = data.groupby("Volume (mL)").mean(numeric_only=True).reset_index()
    volume_data = grouped["Volume (mL)"].values
    potential_data = grouped["Potential (mV)"].values

    # Step 1: Calculate negative value of the derivative of the potential
    dy_dx = -((np.nan_to_num(np.gradient(potential_data, volume_data), nan=0)))

    # Step 2: Define the Cauchy (Lorentzian) function
    def cauchy(x, A, x0, gamma):
        return A / (1 + ((x - x0) / gamma) ** 2)

    window = 0.05  # ±0.05 mL around peak

    # Create mask to keep points outside window OR above a threshold
    x0_guess = volume_data[np.argmax(dy_dx)]
    #threshold = 0.01 * max(dy_dx)  # Keep values above 1% of peak
    #mask = ~((np.abs(volume_data - x0_guess) <= window) & (dy_dx < threshold))

    # Apply mask to both x and y
    #volume_data = volume_data[mask]
    #dy_dx = dy_dx[mask]

    # Step 3: Fit the Cauchy function to the absolute value of the derivative
    initial_guess = [max(dy_dx), x0_guess, 0.01]

    try:
    
        params, covariance = curve_fit(
            cauchy,
            volume_data,
            dy_dx,
            p0=initial_guess)
        

        A_fitted, x0_fitted, gamma_fitted = params
        perr = 2*np.sqrt(np.diag(covariance))  # Parameter standard deviations
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
        ax[1].axvline(x0_fitted, color='green', linestyle=':', label=f"Endpoint (x₀) = {x0_fitted:.4f} ± {2*x0_err:.4f} mL")

    # Print fit parameters to console
    print("\nFit Parameters:")
    print(f"A      = {A_fitted:.2f} ± {A_err:.2f}")
    print(f"x₀     = {x0_fitted:.4f} ± {x0_err:.4f} mL")
    print(f"γ      = {gamma_fitted:.4f} ± {gamma_err:.4f} mL")


    ax[1].set_xlabel("Titrant Volume (mL)")
    ax[1].set_ylabel(r"Absolute Value of Derivative")
    ax[1].set_title("Lorentzian Fit to Derivative of Titration Data")
    ax[1].legend()
    ax[1].grid(True)

    plt.tight_layout()
    plt.show()

    # Ask the user if they'd like to save the plot
    save_plot = input("Would you like to save the plot as a PNG? (y/n): ").strip().lower()
    if save_plot == 'y':
        filename = filepath.split("/")[-1].replace(".csv", "_fit.png")
        try:
            fig.savefig(filename, dpi=300)
            print(f"Plot saved as {filename}")
        except Exception as e:
            print(f"Failed to save figure: {e}")
