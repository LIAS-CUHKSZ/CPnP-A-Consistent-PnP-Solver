# CPnP-A-Consistent-PnP-Solver
The cpp implementation of the CPnP solver is available at https://github.com/bereze/CPnP-cpp-version.

This is the description of the use of the CPnP solver. Paper information:

<pre>
@inproceedings{zeng2023cpnp,
  title={CPnP: Consistent Pose Estimator for Perspective-n-Point Problem with Bias Elimination},
  author={Zeng, Guangyang and Chen, Shiyu and Mu, Biqiang and Shi, Guodong and Wu, Junfeng},
  booktitle={IEEE International Conference on Robotics and Automation (ICRA)},
  pages={},
  year={2023},
  organization={IEEE}
}
</pre>

One can directly run main.m to test the CPnP solver. The main.m mainly includes three parts:
* Generate input data: generate the intrinsics of the camera, the pose of the camera, the coordinates of 3D points, and the noisy 2D projections. This part uses the codes from Copyright (C) <2007>  <Francesc Moreno-Noguer, Vincent Lepetit, Pascal Fua>. Reference paper: Moreno-Noguer F, Lepetit V, Fua P. Accurate Non-Iterative 𝑂 (𝑛) Solution to the P𝑛P Problem[C]//2007 IEEE 11th International Conference on Computer Vision. IEEE, 2007: 1-8.
* Estimate the pose of the camera: estimate the camera pose using CPnP function and calculate the rmse
* Calculate the CRB: calculate the theoretical lower bound for the rmse of any unbiased PnP solver

The inputs and outputs of the CPnP function are as follows:
* Inputs

s - a 3×n matrix whose i-th column is the coordinates (in the world frame) of the i-th 3D point

Psens_2D - a 2×n matrix whose i-th column is the coordinates of the 2D projection of the i-th 3D point
        
fx, fy, u0, v0 - intrinsics of the camera, corresponding to the intrinsic matrix K=[fx 0 u0;0 fy v0;0 0 1]

> **Remark**: the units of Psens_2D, fx, fy, u0, and v0 should be consistent, e.g., all in m, or all in pixels.

* Outputs

R - the estimate of the rotation matrix in the first step
         
t - the estimate of the translation vector in the first step

R_GN - the refined estimate of the rotation matrix with Gauss-Newton iterations

t_GN - the refined estimate of the translation vector with Gauss-Newton iterations
