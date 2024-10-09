# BiSCAN
This is Bidimensional Scan Analyzer (BiSCAN), written in Python! It simultaneously plots the energy profile of a bidimensional relaxed scan and the converged geometries of your system. Anytime you click on a point of the energy profile, the geometry representation updates to show the conformation associated to that point! You only need to provide an output file of a successful bidimensional scan. The program automatically recognizes Orca (version 5) and Gaussian (version 16) output files. It may work with other versions, although it is not guaranteed. 

Installation:

1) You can download the setup.py file and the BiSCAN_module directory and then execute 'pip install .' in the directory where you have both items. That should create the module in your Python installation, which allows you to simply use the code in your Linux terminal as:

biscan output_file.log

where 'output_file.log' is your bidimensional output file. I provide a test file with a bidimensional scan of a water molecule. 




https://github.com/user-attachments/assets/3163d55e-a71e-4fef-a89f-38ef0662ea96




