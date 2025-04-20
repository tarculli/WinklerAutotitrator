'''
This is the script that runs the manual titration mode for Masha's Auto-Titrator. 
To run on Command Prompt, navigate to the "Winkler" folder containing the script, and enter "python manual_titration.py"

Last Update: 04/20/2025
'''

# Importing necessary packages

import winkler_functions as wf
import numpy as np # type: ignore
import time
import matplotlib.pyplot as plt # type: ignore
import re
import csv

# Welcome Message

print(f"""
      
------------------------------------------------------------------------------------------------------------------
\033[1mWelcome to Masha's Auto-Titrator!\033[0m

You are currently running the \033[1mMANUAL TITRATION\033[0m script. 
This mode requires the user to specify the location of the end-point they intend to capture.
Refer to the instrument protocol to ensure that it has been properly set-up.
------------------------------------------------------------------------------------------------------------------
      
      """)
time.sleep(3.5)
input(f"Press \033[1;32mENTER\033[0m when you are ready to begin the titration process...\n")

name = input(f"\033[1mWho is operating the auto-titrator?\033[0m\n")
sample_number = input(f"\n\033[1mWhat is the sample ID? (ex. cruise_bottlenumber)\033[0m\n")
syringe_D = float(input(f"\n\033[1mWhat is the loaded syringe diameter (mm)?\033[0m\n"))
ramp_phase = float(input(f"\n\033[1mHow long do you want to ramp towards the capture interval?\033[0m\n"))
capture_phase = float(input(f"\n\033[1mHow broad do you want to set the end-point capture interval?\033[0m\n"))

# Writing Pump Program

program_string = f"""
DIA {syringe_D}

PHN 1
FUN RAT
RAT 0.5 MM
VOL {ramp_phase}
DIR INF

PHN 2
FUN RAT
RAT 0.02 MM
VOL {capture_phase}
DIR INF

PHN 3
FUN RAT
RAT 0.5 MM
VOL {ramp_phase}
DIR INF

PHN4
FUN BEP

PHN 5
FUN STP
"""

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

# Sending Program to Pump

for line in program_string.split('\n'):
    line = line.strip()  # remove leading/trailing whitespaces
    if line:  # skip empty lines
        wf.send_pump_command(line, pump) # send the command line to the pump

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

volume_data = [] 
mV_data = []

try: 
    while True:
        wf.send_pump_command("DIS", pump)
        volume = re.search(r"I(\d{1,2}\.\d{3})W",pump.readline().decode().strip()).group(1)
        meter_data = meter.readline().decode("utf-8").strip()

        if float(volume) >= (2*ramp_phase+capture_phase): 
            print("\n\033[1;32mPump has finished its program. Stopping data collection.\033[0m")
            break  # Exit the loop

        if meter_data:
            mV = float(re.search(r"mV,(-?\d+\.\d+),mV", meter_data).group(1))
            print(f"Pump Output: {volume} mL \nMeter Output: {mV} mV")
            volume_data.append(float(volume))
            mV_data.append(mV)
            plt.scatter(float(volume), mV, color='b')
            plt.draw()
            plt.pause(0.1)
        time.sleep(0.5)

# An option (CTRL+C) to kill the data collection loop if something goes wrong or data has been collected to a satifactory degree

except KeyboardInterrupt:
    wf.send_pump_command("STP", pump)
    wf.close_connections(meter, pump)

# Save data to CSV

to_save = input(f"\n\033[1mPress enter to save the data from this run...").strip().lower()


# Metadata
metadata = [
f"# Titration Data",
f"# Sample ID: {sample_number}",
f"# Titration Date: {time.strftime('%Y-%m-%d %H:%M:%S')}",
f"# Operator: {name}",
f"# Notes: Titration using manual mode.",
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
