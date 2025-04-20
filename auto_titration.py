'''
This is the script that runs the automatic titration mode for Masha's Auto-Titrator. 
To run on Command Prompt, navigate to the "Winkler" folder containing the script, and enter "python auto_titration.py"

Last Update: 04/20/2025
'''

# Importing necessary packages

import winkler_functions as wf
import numpy as np # type: ignore
import time
import matplotlib.pyplot as plt # type: ignore
import re
import csv
import warnings
warnings.simplefilter("ignore", category=RuntimeWarning)

# Welcome Message

print(f"""
      
------------------------------------------------------------------------------------------------------------------
\033[1mWelcome to Masha's Auto-Titrator!\033[0m

You are currently running the \033[1mAUTOMATIC TITRATION\033[0m script. 
This mode dispenses titrant at a rate that is inversely proportional to the magnitude of change in potential.
Refer to the instrument protocol to ensure that it has been properly set-up.
------------------------------------------------------------------------------------------------------------------
      
      """)
time.sleep(3.5)
input(f"\nPress \033[1;32mENTER\033[0m when you are ready to begin the titration process...\n")

name = input(f"\033[1mWho is operating the auto-titrator?\033[0m\n")
sample_number = input(f"\n\033[1mWhat is the sample ID? (ex. cruise_bottlenumber)\033[0m\n")
syringe_D = float(input(f"\n\033[1mWhat is the loaded syringe diameter (mm)?\033[0m\n"))
end_volume = float(input(f"\n\033[1mAt what volume do you want to terminate the titration process (mL)?\033[0m\n"))


# Opening Serial Connections to Pump and Meter

input(f"""
-----------------------------------------------------------
Press \033[1;32mENTER\033[0m to open serial connections...
-----------------------------------------------------------
""")

pump = wf.connect_pump_MIMSy()
meter = wf.connect_meter_MIMSy()

wf.send_pump_command("RESET", pump)

input(f"\n\033[32mConnections active...\033[0m\n\nPress \033[1;32mENTER\033[0m to continue.\n")

wf.send_pump_command(f"DIA {syringe_D}", pump)
wf.send_pump_command(f"RAT 0.2 MM", pump)
wf.send_pump_command("STP", pump)

print("Titration setup complete.")
time.sleep(2.0)
print("""

-------------------------------------------------------------------
IMPORTANT:
Before moving on, clear the data log on your Orion Star A321 meter.
-------------------------------------------------------------------

Follow these steps:

    1) Select \"log view\" by pressing f3.
    2) Use the select option (f2) to navigate to the \"Data Log\" menu.
    3) Open the \"options\" menu (f3).
    4) Select \"Log Clear\" and press f2 to confirm.
    5) Press f3 to begin entering the password (111111) and f2 to submit the password.
    6) Continue to the next step.

""")
time.sleep(2.0)
input("\n\n\nWhen you are ready, press \"measure\" on the meter, then quickly press \033[1;32mENTER\033[0m to continue...")


# Titration and Data Collection

sensitivity = 25
volume_data = [] 
mV_data = []
thresholds = sorted([sensitivity, sensitivity*2.5, sensitivity*3.5, sensitivity*5], reverse=True)


try: 
    while True:

        # Looping to query volume dispensed from the pump and simultaneous potential readout from the meter 

        wf.send_pump_command("DIS", pump)
        volume = re.search(r"I(\d{1,2}\.\d{3})W",pump.readline().decode().strip()).group(1)
        meter_data = meter.readline().decode("utf-8").strip()

        # Conditionally ending the data collection loop once enough volume has been dispensed

        if float(volume) >= (end_volume): 
            wf.send_pump_command("STP", pump)
            wf.send_pump_command("BEP", pump)
            print("\n\033[1;32mPump has finished its program. Stopping data collection.\033[0m")
            break
    
        # Storing the data if it exists and calculating the derivatives to inform pumping rate

        if meter_data:
            mV = float(re.search(r"mV,(-?\d+\.\d+),mV", meter_data).group(1))
            if volume_data:
                x = [volume_data[-1],float(volume)]
            else:
                x = [0, float(volume)]
            if mV_data:
                y = [mV_data[-1],mV]
            else:
                y = [0, float(volume)]
            dx = abs(np.diff(x))
            if dx < 0.01:
                dx = 0.01
            dy = abs(np.diff(y))
            dydx = dy/dx
            print(f"Pump Output: {volume} ml \nMeter Output: {mV} mV \nDifference Parameter: {dydx}")
            volume_data.append(float(volume))
            mV_data.append(mV)

            # Conditionally altering the pumping rate based on the current derivative of the output data

            if dydx >= thresholds[0]:
                wf.send_pump_command("STP", pump)
                wf.send_pump_command(f"RAT 0.02 MM", pump)
                wf.send_pump_command("RUN", pump)
            elif dydx >= thresholds[1]:
                wf.send_pump_command("STP", pump)
                wf.send_pump_command(f"RAT 0.05 MM", pump)
                wf.send_pump_command("RUN", pump)
            elif dydx >= thresholds[2]:
                wf.send_pump_command("STP", pump)
                wf.send_pump_command(f"RAT 0.1 MM", pump)
                wf.send_pump_command("RUN", pump)
            elif dydx >= thresholds[3]:
                wf.send_pump_command("STP", pump)
                wf.send_pump_command(f"RAT 0.2 MM", pump)
                wf.send_pump_command("RUN", pump)
            else:
                wf.send_pump_command("STP", pump)
                wf.send_pump_command(f"RAT 0.5 MM", pump)
                wf.send_pump_command("RUN", pump)
            
            # Updating the live plot with the latest data point

            plt.scatter(float(volume), mV, color='b', s=8)
            plt.draw()
            plt.pause(0.1)
            
        time.sleep(0.5)

# An option (CTRL+C) to kill the data collection loop if something goes wrong or data has been collected to a satifactory degree

except KeyboardInterrupt:
    wf.send_pump_command("STP", pump)
    wf.close_connections(meter, pump)

# Save data to CSV

to_save = input(f"Press enter to save the data from this run...").strip().lower()


# Metadata
metadata = [
f"# Titration Data",
f"# Sample ID: {sample_number}",
f"# Titration Date: {time.strftime('%Y-%m-%d %H:%M:%S')}",
f"# Operator: {name}",
f"# Notes: Titration using automatic mode.",
f"# Columns: Volume (mL), Potential (mV)"]

filename = f"Data/{sample_number}_titration_{time.strftime('%Y%m%d')}.csv"

with open(filename, mode='w', newline='') as file:
    
    # Write metadata as comments
    for line in metadata:
        file.write(line + "\n")  # Writing plain text comments

    writer = csv.writer(file)
    writer.writerow(["Volume (mL)", "Potential (mV)"])  # Column headers
    writer.writerows(zip(volume_data, mV_data))  # Write data rows

print(f"\033[1mData saved to {filename}\033[0m")
