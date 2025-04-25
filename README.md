# WinklerAutotitrator
A repository for the functions and scripts needed to run an automatic titration on Masha's Auto-Titrator.

These are the Python scripts necessary to run titrations and process titration data with the titration system built by Toby Arculli PO'25 in the Prokopenko Lab.

The "auto_titration.py" script is used for automatic titrations and relies on the serial communication between the pump, meter, and control computer to cycle between pumping rates based on the immediate "derivative" of the titration data. To run this script, navigate to the Winkler folder in the MIMSy laptop, run "python auto_titration.py", and follow the prompts. If you need to set the system up with a new computer, make sure to alter the commands that connect the pump and meter so that they reflect the correct RS232 ports (see "winkler_functions.py").

The "manual_titration.py" script runs similarly, but it needs user input to define the regions whee it pumps rapidly or slowly.

Once data is collected and stored, you can run the "lorentzian_end_point.py" script to plot the data and extract and end-point fit from a Lorentzian model on the first derivatives to the titration data. To do so, run "python lorentzian_end+point.py" and enter the filepath for the data you wish to plot when prompted.

"winkler_functions.py" is a collection of some useful and some archived functions that were developed for the autotitrator when we were building the device.
