# ttk4900_master_thesis
This code finds a controller for the Butterfly Robot. The periodic Riccati differential equation is solved by the function sdp_riccati in the folder sdp.This code needs two dependencies YALMIP and SDPT3.

In the class butterfly_robot all the calculations for the butterfly robot
are done. Instantiate this with the command bf = butterfly_robt(true/false,true/false);
It takes two options, the first determines whether the LQR is calculated, the second determines whether an integral term should be added to the LQR(not working well at the moment.). 
The LQR can be tuned in the private variables for the class.
If an instance of the butterfly_robot exists with the name bf, the script generate_data can be run to create spline interpolations of the Riccati solution, the nominal trajectory and the shape of the frame. These are needed for the simulation of the system found in butterfly_robot_matlab_to_c.slx.These are also needed for the function controller.m.
To generate the C++ code needed to run the real Rutterfly Robot use matlab coder to generate C++ code from the matlab function controller.m
