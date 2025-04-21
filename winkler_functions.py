import numpy as np # type:ignore
import pandas as pd # type:ignore
import matplotlib.pyplot as plt # type:ignore
import serial # type:ignore
import time

'''

Winkler Titration Model

Description: This function yields a modeled titration based on the Nernst Equations for the Winkler Method

Inputs: A theoretical dissolved oxygen concentration in micromolar units (o2), the molarity of the thiosulfate titrant (M_thio),
the volume of the sample in liters (V_sample), the temperature at which the theoretical reaction occurs in kelvin (T), 
and a boolean argument for the output of a corresponding plot (if desired, plotting=True)

Outputs: The equivalence volume (ie. endpoint volume) of the modeled titration (V_equivalence)in L, a volume array in mL (volume),
and a potential array in Volts (potential)

To-do: Make units uniform and more intuitive!

Last Update: 2/27/25 -- Toby Arculli

Example Usage: 
v_eq, v, p = winkler_model(301,0.2,0.1,298) 

'''

def winkler_model(o2, M_thio, V_sample, T, plotting=False):

  # CONSTANTS
  R = 1.987*10**-3 # gas constant kcal/deg*mol
  F = 23.03 # kcal/V
  E_0_thio = 198*10**-6 # standard potential in volts
  n_thio = 2 # number of electrons transferred in moles
  E_0_io = 0.53 # standard potential in volts
  n_io = 2 # number of electrons transferred in moles

  # MODEL INPUTS

  T = T # temperature in K
  V_sample = V_sample # volume of sample in liters
  M_thiosulfate_titrant = M_thio # molarity of thiosulfate
  O2_concentration = o2*10**-6 # molarity of oxygen

  # EQUIVALENCE VOLUME
  # calculating the volume of titrant needed to reach the equivalence point
  V_equivalence = (2*O2_concentration*V_sample)*2/(M_thiosulfate_titrant)

  # TITRATION CUMULATIVE VOLUME
  V_titrant_added = np.linspace(0,2*V_equivalence,5000) # Creating a range of values for the amount of titrant added

  # IODINE-IODIDE CONCENTRATIONS
  # moles of iodine with the addition of titrant (essentially intial iodine - a half mole for every mole of thiosulfate added)
  moles_iodine = 2*O2_concentration*V_sample-0.5*M_thiosulfate_titrant*V_titrant_added
  # moles of iodide with the addition of titrant (stoichiometrically twice the amount of iodine consumed in the reaction)
  moles_iodide = 2*(2*O2_concentration*V_sample-moles_iodine)

  M_iodine = moles_iodine/(V_sample+V_titrant_added) # concentration of iodine
  M_iodide = moles_iodide/(V_sample+V_titrant_added) # concentration of iodide

  Q_io = M_iodide**2/M_iodine # ratio of molarities

  # TETRATHIONIDE-THIOSULFATE CONCENTRATIONS
  # moles of tetrathionide (equivalent to the moles of reacted iodine)
  moles_tetrathionide = 2*O2_concentration*V_sample-moles_iodine
  # moles of thiosulfate (zero until all iodine is reacted, then equal to the product of added volume and molarity)
  moles_thiosulfate = M_thiosulfate_titrant*V_titrant_added-2*2*O2_concentration*V_sample
  M_tetrathionide = moles_tetrathionide/(V_sample+V_titrant_added) # concentration of tetrathionide
  M_thiosulfate = moles_thiosulfate/(V_sample+V_titrant_added) # concentration of thiosulfate

  Q_thio = M_thiosulfate**2/M_tetrathionide # ratio of molarities

  # NERNST EQUATIONS
  E_thio = np.where(Q_thio > 0, E_0_thio - ((R*T)/(n_thio*F))*np.log(Q_thio), 0)
  E_io = np.where(Q_io > 0, E_0_io - ((R*T)/(n_io*F))*np.log(Q_io), 0)

  E_equivalence_a = (E_0_thio/2 + E_0_io/2 - ((R*T)/(4*F))*np.log(8*M_iodine*M_iodide))
  upper = E_equivalence_a[np.where(V_titrant_added < V_equivalence)[0][-1]]
  E_equivalence_b = (E_0_thio/2 + E_0_io/2 - ((R*T)/(4*F))*np.log(4*M_thiosulfate*M_iodide))
  lower = E_equivalence_b[np.where(V_titrant_added > V_equivalence)[0][1]]
  E_equivalence = (upper+lower)/2

  iodine_mask = V_titrant_added < V_equivalence
  thiosulfate_mask = V_titrant_added >= V_equivalence

  volume = np.concatenate((V_titrant_added[iodine_mask] * 10**3,
                                 V_titrant_added[thiosulfate_mask] * 10**3))
  potential = np.concatenate((E_io[iodine_mask],
                                 E_thio[thiosulfate_mask]))



  if plotting:

    # PLOTTING RELEVANT REGIONS
    plot = plt.figure()
    plt.scatter(volume,potential,label='Predicted Potential', s=1)
    plt.scatter(V_equivalence*(10**3),E_equivalence,label='Equivalence Point',color='red', marker='X')
    plt.vlines(x=V_equivalence*(10**3), ymin=0, ymax=E_equivalence, color='red',linestyle='--',label='End-Point Volume', zorder=2)
    plt.title(f'Model of Potential Response to the Addition of Thiosulfate to a {V_sample} L sample\n'
    + r'$(O_2 = ' + str(o2) + r' \, \mu mol/L, \, M_{thio} = ' + str(M_thio) + r' \, M, \, T = ' + str(T) + r' \, K)$')
    plt.ylabel('Voltage (V)')
    plt.xlabel('Volume of titrant added (mL)')
    plt.ylim(bottom=0)
    plt.legend()
    plt.show()

    return V_equivalence, volume, potential

  else:
    return V_equivalence, volume, potential




    
'''

Derivative-based End-point Determination

Description: This function calculates the end-point based on the largest dericative value in the titration data

Inputs: Volume data (volume_data) in mL and potential data (potential_data) in mV

Outputs: Prints the volume associated with the end-point

To-do: N/A

Last Update: 2/27/25 -- Toby Arculli

Example Usage: N/A -- quite self-explanatory... (ARCHIVED)

'''

def find_endpoint(volume_data, potential_data):
    dy_dx = np.gradient(potential_data, volume_data)  # Compute first derivative
    min_index = np.argmin(dy_dx) # Find index of maximum derivative
    max_derivative = dy_dx[min_index]  # Max derivative value
    corresponding_volume = volume_data[min_index]  # Corresponding volume

    print(f"Equivalence Volume: {corresponding_volume}")




'''

A Simple Derivative Calculator

Description: This function calculates the derivatives of given x and y data and yields a plot of the derivative vs. x-data
if plotting=True

Inputs: X data, Y data, plotting option

Outputs: Derivatives (dy_dx) and an optional plot of dy_dx vs. x

To-do: N/A

Last Update: 2/27/25 -- Toby Arculli

Example Usage: Again... quite straightforward

'''

def find_gradient(volume_data, potential_data, plotting=False):
    dy_dx = np.gradient(potential_data, volume_data)  # Compute first derivative
    if plotting:
      plt.plot(volume_data,dy_dx)
    else:
      return dy_dx




'''

Ideal Pump-rate Optimization

Description: Adjusts the pumping rate based on the steepness of the derivative and calculate total pumping time.
Inputs:

  volume_data (array): Volume of titrant added (mL)
  potential_data (array): Corresponding voltage (V)
  R_max (float): Maximum pump rate (microliters/min)
  R_min (float): Minimum pump rate (mL/min) to prevent stopping
  k (float): Scaling factor to control sensitivity
  plotting (bool): Whether to plot the adjusted pump rate

Outputs:

  pump_rates (array): Adjusted pump rates for each volume step
  total_pumping_time (float): Total time required for titration (minutes)

To-do: N/A

Last Update: 2/27/25 -- Toby Arculli 

Example Usage: (ARCHIVED)

'''

def optimize_pump_rate(volume_data, potential_data, R_max=1, R_min=0.003, k=0.1, plotting=False):
  
    dy_dx = np.gradient(potential_data, volume_data)  # Compute first derivative
    
    # Adjust pump rate based on the derivative steepness
    pump_rates = R_max / (1 + k * np.abs(dy_dx))  # Compute initial rates

    # Ensure pump rate does not fall below R_min
    pump_rates = np.maximum(pump_rates, R_min)

    # Compute volume increments (∆V) in mL
    dV = np.diff(volume_data)

    # Compute time increments (∆t = ∆V / Pump Rate)
    dt = dV / pump_rates[:-1]  # Use pump rate at corresponding step

    # Total pumping time (sum of all ∆t)
    total_pumping_time = np.sum(dt)

    if plotting:
        fig, ax1 = plt.subplots()

        # Plot potential vs volume
        ax1.set_xlabel("Volume of Titrant (mL)")
        ax1.set_ylabel("Potential (V)", color="blue")
        ax1.scatter(volume_data, potential_data, label="Potential", color="blue", s=10)
        ax1.tick_params(axis="y", labelcolor="blue")

        # Plot adjusted pump rate on a second y-axis
        ax2 = ax1.twinx()
        ax2.set_ylabel("Pump Rate (mL/min)", color="red")
        ax2.plot(volume_data, pump_rates, label="Pump Rate", color="red", linestyle="dashed")
        ax2.tick_params(axis="y", labelcolor="red")

        plt.title(f"Pump Rate Optimization (Total Time: {total_pumping_time:.2f} min)")
        plt.legend()
        plt.show()

    return pump_rates, total_pumping_time




'''
Simplified Pump Rates for Easy Programming
(ARCHIVED)
'''

def realistic_pump_rates(pump_rates):
  pump_rates[0:1400] = max(pump_rates)
  for i in range(len(pump_rates)):
    if pump_rates[i] > 0.999*max(pump_rates):
      pump_rates[i] = max(pump_rates)
    else:
      pump_rates[i] = 0.02

  new_pump_rates = pump_rates

  step_indices = np.where(np.diff(new_pump_rates) != 0)[0]

  return new_pump_rates, step_indices


def get_transition_volumes(pump_rates, volume_data):
    """Finds volumes where pump rates change."""
    new_pump_rates, step_indices = realistic_pump_rates(pump_rates)

    transition_volumes = [volume_data[i] for i in step_indices]

    return transition_volumes


'''
A function to write pumping programs given a Winkler model, ideal pumping rates, and syringe diameter.
(ARCHIVED)
'''

def winkler_pumping_program(pump_rates, syringe_diameter, volume_data, step_indices):
    step_indices = np.array(step_indices).astype(int)  # Ensure indices are integers

    transition_volumes = volume_data[step_indices]
    rates = pump_rates[step_indices]
    #{(transition_volumes[1] - transition_volumes[0]):.3f}
    program_string = f"""
    DIA {syringe_diameter}

    PHN 1
    FUN RAT
    RAT {rates[0]:.2f} MM
    VOL {transition_volumes[0]:.3f}
    DIR INF

    PHN 2
    FUN RAT
    RAT {rates[1]:.2f} MM
    VOL 0.15
    DIR INF

    PHN 3
    FUN RAT
    RAT {rates[0]:.2f} MM
    VOL {transition_volumes[0]:.3f}
    DIR INF

    PHN4
    FUN BEP

    PHN 5
    FUN STP
    """

    return program_string






'''
Establish Serial Connection to OrionStar A321 meter

Description:

Inputs: None

Outputs: None

To-do: N/A

Last Update: 04/15/25

Example Usage: Just run the command once physical connections are active.
'''

def connect_meter():
  import serial.tools.list_ports # type:ignore
    
  # List available COM ports
  ports = [port.device for port in serial.tools.list_ports.comports()]
  print("Available COM ports:\n", ports)

  port = input("Enter the port name for the pump connection (see protocol for Mac vs. Windows naming): ")
  port_name = f'{port}'

  try:
      meter = serial.Serial(port_name, 19200, parity=serial.PARITY_NONE, bytesize=8, stopbits=1, timeout=0.1, xonxoff=True, rtscts=False)
      print(f"Connected successfully to {port_name}")
      return meter
  except serial.SerialException as e:
      print(f"Failed to connect to {port_name}: {e}")
      return None
  
# Here is a version specific to the MIMSy computer and 2025 protocol
  
def connect_meter_MIMSy():

  port_name = 'COM3'
  try:
      pump = serial.Serial(port_name, 19200, parity=serial.PARITY_NONE, bytesize=8, stopbits=1, timeout=0.1, xonxoff=False, rtscts=False)
      print(f"Connected successfully to {port_name}")
      return pump
  except serial.SerialException as e:
      print(f"Failed to connect to {port_name}: {e}")
      return None



'''
Establish Serial Connection to NE-1000 pump

Description:

Inputs: None

Outputs: None

To-do: N/A

Last Update: 04/15/25

Example Usage: Just run the command once physical connections are active.
'''

def connect_pump():
  import serial.tools.list_ports # type:ignore
    
  # List available COM ports
  ports = [port.device for port in serial.tools.list_ports.comports()]
  print("Available COM ports:\n", ports)

  port = input("Enter the port name for the pump connection (see protocol for Mac vs. Windows naming): ")
  port_name = f'{port}'

  try:
      pump = serial.Serial(port_name, 19200, parity=serial.PARITY_NONE, bytesize=8, stopbits=1, timeout=0.1, xonxoff=False, rtscts=False)
      print(f"Connected successfully to {port_name}")
      return pump
  except serial.SerialException as e:
      print(f"Failed to connect to {port_name}: {e}")
      return None
  
# Here is a version of the command specific to the MIMSy laptop and 2025 protocol
  
def connect_pump_MIMSy():

  port_name = 'COM6'
  try:
      pump = serial.Serial(port_name, 19200, parity=serial.PARITY_NONE, bytesize=8, stopbits=1, timeout=0.1, xonxoff=False, rtscts=False)
      print(f"Connected successfully to {port_name}")
      return pump
  except serial.SerialException as e:
      print(f"Failed to connect to {port_name}: {e}")
      return None

'''
Sending a command to the pump
'''

def send_pump_command(command, pump_serial_object):
    pump_serial_object.write((command + '\r\n').encode())  # send commands as encoded strings # type: ignore
    time.sleep(0.1)  # brief delay allowing time for the pump to process the command




'''
Close all connections
'''

def close_connections(meter_serial_object, pump_serial_object):
   meter_serial_object.close()   # type: ignore
   pump_serial_object.close()    # type: ignore
   print("Closing serial connection to pump and meter...")

  

'''
Lorentzian Fit to find end-point

Description:
This function compites the end-point by fitting a lorentzian function to the first derivative of
the potentiometric titration data.

Inputs: None

Outputs: None

To-do: N/A

Last Update: 04/15/25

Example Usage: Just run the command once physical connections are active.
'''

def lorentzian_fit(data_file):
    # Read the CSV file while skipping comment lines
    try:
        data = pd.read_csv(data_file, comment='#')  # Skip metadata lines
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

        # Step 2: Define the Cauchy function
        def cauchy(x, A, x0, gamma):
            return A / (1 + ((x - x0) / gamma) ** 2)

        # Step 3: Fit the Cauchy function to the absolute value of the derivative
        initial_guess = [max(dy_dx), volume_data[np.argmax(dy_dx)], 0.001]

        try:
            params, covariance = curve_fit(cauchy, volume_data, dy_dx, p0=initial_guess) # type:ignore
            A_fitted, x0_fitted, gamma_fitted = params
            x0_err = np.sqrt(np.diag(covariance))[1]  # Standard deviation of the fitted x0
            fitted_curve = cauchy(np.linspace(0, max(volume_data), 10000), *params)
        except Exception as e:
            print(f"Curve fitting failed: {e}")
            params = None

        # Create subplots
        fig, ax = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

        # Top Panel: Raw Titration Data
        ax[0].scatter(volume_data, potential_data, color='blue', label="Raw Data", s=12)
        ax[0].set_ylabel("Potential (mV)")
        ax[0].set_title(f"{data_file[:-4]}")

        # Add a line for the endpoint location
        if params is not None:
            ax[0].axvline(x0_fitted, color='green', linestyle=':', label=f"Endpoint (x0) = {x0_fitted:.3f} mL")

        ax[0].legend()
        ax[0].grid(True)

        # Bottom Panel: Derivative and Cauchy Fit
        ax[1].scatter(volume_data, dy_dx, label="First Derivative", color='blue', s=12)

        if params is not None:
            ax[1].plot(np.linspace(0, max(volume_data), 10000), fitted_curve, label="Lorentzian Fit", color='red', linestyle='--')
            ax[1].axvline(x0_fitted, color='green', linestyle=':', label=f"Endpoint (x0) = {x0_fitted:.3f} ± {2*x0_err:.5f} mL")

        ax[1].set_xlabel("Titrant Volume (mL)")
        ax[1].set_ylabel(r"Absolute Value of Derivative")
        ax[1].set_title("Lorentzian Fit to Derivative of Titration Data")
        ax[1].legend()
        ax[1].grid(True)

        # Ensure the figure is rendered before accessing the buffer
        fig.canvas.draw()

        return fig  # Return the figure object
