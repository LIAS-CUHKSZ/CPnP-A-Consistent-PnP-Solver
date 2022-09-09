# CPnP_A-Consistent-PnP-Solver
The is the description of the use of the CPnP solver. 

One can directly run main.m to test the CPnP solver. The main.m mainly includes three parts:
* Generate input data: generate the intrinsics of the camera, the pose of the camera, the coordinates of 3D points, and the noisy 2D projections. This part uses the codes from Copyright (C) <2007>  <Francesc Moreno-Noguer, Vincent Lepetit, Pascal Fua>. Paper title: Accurate Non-Iterative 𝑂 (𝑛) Solution to the P𝑛P Problem.
* Estimate the pose of the camera: estimate the camera pose using CPnP function and calculate the rmse
* Calculate the CRB: calculate the theoretical lower bound for the rmse of any unbiased PnP solver
