# CPnP: a consistent PnP solver
# Inputs: s - a 3×n matrix whose i-th column is the coordinates (in the world frame) of the i-th 3D point
# Psens_2D - a 2×n matrix whose i-th column is the coordinates of the 2D projection of the i-th 3D point
# fx, fy, u0, v0 - intrinsics of the camera, corresponding to the intrinsic matrix K=[fx 0 u0;0 fy v0;0 0 1]

# Outputs: R - the estimate of the rotation matrix in the first step
# t - the estimate of the translation vector in the first step
# R_GN - the refined estimate of the rotation matrix with Gauss-Newton iterations
# t_GN - the refined estimate of the translation vector with Gauss-Newton iterations
# Copyright <2022> <Guangyang Zeng, Shiyu Chen, Biqiang Mu, Guodong Shi, Junfeng Wu>
# Guangyang Zeng, SLAMLab-CUHKSZ, September 2022
# zengguangyang@cuhk.edu.cn, https://github.com/SLAMLab-CUHKSZ 
# paper link: https://arxiv.org/abs/2209.05824



import numpy as np
from numpy import linalg
from scipy.linalg import expm, eigh, eig, svd

def CPnP(s, Psens_2D, fx, fy, u0, v0):
    N = s.shape[1]
    bar_s = np.mean(s, axis=1).reshape(3,1)
    Psens_2D = Psens_2D - np.array([[u0],[v0]])
    obs = Psens_2D.reshape((-1, 1), order="F")
    pesi = np.zeros((2*N,11))
    G = np.ones((2*N,1))
    W = np.diag([fx, fy])
    M = np.hstack([np.kron(bar_s.T, np.ones((2*N,1))) - np.kron(s.T, np.ones((2,1))), np.zeros((2*N, 8))])

    for k in range(N):
        pesi[[2*k],:] = np.hstack([-(s[0,k]-bar_s[0]) * obs[2*k], -(s[1,k]-bar_s[1]) * obs[2*k], -(s[2,k]-bar_s[2]) * obs[2*k], (fx * s[:,[k]]).T.tolist()[0], fx, 0, 0, 0, 0])
        pesi[[2*k+1],:] = np.hstack([-(s[0,k]-bar_s[0]) * obs[2*k+1], -(s[1,k]-bar_s[1]) * obs[2*k+1], -(s[2,k]-bar_s[2]) * obs[2*k+1], 0, 0, 0, 0, (fy * s[:,[k]]).T.tolist()[0], fy])

    J = np.dot(np.vstack([pesi.T, obs.T]), np.hstack([pesi, obs])) / (2*N)
    delta= np.vstack([np.hstack( [np.dot(M.T, M), np.dot(M.T, G)]), np.hstack( [np.dot(G.T,M), np.dot(G.T, G)] )]) / (2*N)

    w, D = eig(J, delta)    
    sigma_est = min(abs(w))

    est_bias_eli= np.dot(np.linalg.inv((np.dot(pesi.T,pesi)- sigma_est*(np.dot(M.T, M)))/(2*N)) , (np.dot(pesi.T,obs) - sigma_est * np.dot(M.T, G))/ (2*N) )
    bias_eli_rotation = np.vstack([est_bias_eli[3:6].T, est_bias_eli[7:10].T, est_bias_eli[0:3].T])
    bias_eli_t= np.hstack([est_bias_eli[6],  est_bias_eli[10], 1-bar_s[0]*est_bias_eli[0]-bar_s[1]*est_bias_eli[1]-bar_s[2]*est_bias_eli[2]]).T
    normalize_factor= np.linalg.det(bias_eli_rotation) ** (1/3)
    bias_eli_rotation=bias_eli_rotation / normalize_factor
    t = bias_eli_t/normalize_factor

    U, x, V = svd(bias_eli_rotation)
    V = V.T

    RR = np.dot(U, np.diag([1, 1, np.linalg.det(np.dot(U, V.T))]))
    R = np.dot(RR, V.T)

    E = np.array([[1, 0, 0],[0, 1, 0]])
    WE = np.dot(W, E)
    e3 = np.array([[0],[0],[1]])
    J = np.zeros((2*N, 6))

    g = np.dot(WE, np.dot(R, s) + np.tile(t,N).reshape(N,3).T)
    h = np.dot(e3.T, np.dot(R, s) + np.tile(t,N).reshape(N,3).T)

    f = g/h
    f =  f.reshape((-1, 1), order="F")
    I3 = np.diag([1, 1, 1])

    for k in range(N):
        J[[2*k, 2*k+1],:] = np.dot((WE * h[0, k] - g[:,[k]]* e3.T), np.hstack([s[1,k]*R[:,[2]]-s[2,k]*R[:,[1]], s[2,k]*R[:,[0]]-s[0,k]*R[:,[2]], s[0,k]*R[:,[1]]-s[1,k]*R[:,[0]], I3])) / h[0, k]**2

    initial = np.hstack([np.zeros((3)), t.tolist()]).reshape(6,1)
    results = initial + np.dot(np.dot(np.linalg.inv(np.dot(J.T, J)), J.T), (obs - f))
    X_GN = results[0:3]
    t_GN = results[3:6]
    Xhat = np.array([
    [0, -X_GN[2], X_GN[1]], 
    [X_GN[2], 0, -X_GN[0]],
    [-X_GN[1], X_GN[0], 0]
    ])
    R_GN = np.dot(R, expm(Xhat))

    return R,t,R_GN,t_GN



