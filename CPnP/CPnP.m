% CPnP: a consistent PnP solver
% Inputs: s - a 3*n matrix whose i-th column is the coordinates (in the world frame) of the i-th 3D point
%        Psens_2D - a 2*n matrix whose i-th column is the coordinates of the 2D projection of the i-th 3D point
%        fx, fy, u0, v0 - intrinsics of the camera, corresponding to the intrinsic matrix K=[fx 0 u0;0 fy v0;0 0 1]
%
%Outputs: R - the estimate of the rotation matrix in the first step
%         t - the estimate of the translation vector in the first step
%         R_GN - the refined estimate of the rotation matrix with Gauss-Newton iterations
%         t_GN - the refined estimate of the translation vector with Gauss-Newton iterations

% Copyright <2022>  <Guangyang Zeng, Shiyu Chen, Biqiang Mu, Guodong Shi, Junfeng Wu>
% Guangyang Zeng, SLAMLab-CUHKSZ, September 2022
% zengguangyang@cuhk.edu.cn, https://github.com/SLAMLab-CUHKSZ 
% paper information: Guangyang Zeng, Shiyu Chen, Biqiang Mu, Guodong Shi, and Junfeng Wu. 
%                    CPnP: Consistent Pose Estimator for Perspective-n-Point Problem with Bias Elimination, 
%                    IEEE International Conference on Robotics and Automation (ICRA), London, UK, May 2023.

function [R,t,R_GN,t_GN]=CPnP(s,Psens_2D,fx,fy,u0,v0)
N=length(s);
bar_s=sum(s,2)/N; 
Psens_2D=Psens_2D-[u0;v0];
obs=Psens_2D(:);  
pesi=zeros(2*N,11); 
G=ones(2*N,1);
W=diag([fx fy]);
M=[kron(bar_s',ones(2*N,1))-kron(s',ones(2,1)) zeros(2*N,8)];
for k=1:N
    pesi(2*k-1,:)=[-(s(1,k)-bar_s(1))*obs(2*k-1) -(s(2,k)-bar_s(2))*obs(2*k-1) -(s(3,k)-bar_s(3))*obs(2*k-1) fx*s(:,k)' fx 0 0 0 0];
    pesi(2*k,:)=[-(s(1,k)-bar_s(1))*obs(2*k) -(s(2,k)-bar_s(2))*obs(2*k) -(s(3,k)-bar_s(3))*obs(2*k) 0 0 0 0 fy*s(:,k)' fy];
end
%% noise variance estimation
J=[pesi';obs']*[pesi obs]/(2*N);
delta=[M'*M M'*G;G'*M G'*G]/(2*N);
[~,D]=eig(J,delta);
root=diag(abs(D));
sigma_est=min(root);

%% First step: bias-eliminated solution
est_bias_eli=((pesi'*pesi-sigma_est*(M'*M))/(2*N))\((pesi'*obs-sigma_est*M'*G)/(2*N));
bias_eli_rotation=[est_bias_eli(4:6)';est_bias_eli(8:10)';est_bias_eli(1:3)'];
bias_eli_t=[est_bias_eli(7) est_bias_eli(11) 1-bar_s(1)*est_bias_eli(1)-bar_s(2)*est_bias_eli(2)-bar_s(3)*est_bias_eli(3)]';
normalize_factor=(det(bias_eli_rotation))^(1/3);
bias_eli_rotation=bias_eli_rotation/normalize_factor; 
t=bias_eli_t/normalize_factor; 
[U,~,V] = svd(bias_eli_rotation);
R=U*diag([1 1 det(U*V')]')*V';

%% Second step: GN iterations
E=[1 0 0;0 1 0];
WE=W*E;
e3 = [0;0;1];
J = zeros(2*N, 6); 
g=WE* (R*s+t);
h=e3'* (R*s+t);
f=g./h;
f=f(:);
I3=eye(3);
for k = 1:N  
    J(2*k-1:2*k,:) = (((WE * h(k) - g(:,k)* e3') * [s(2,k)*R(:,3)-s(3,k)*R(:,2) s(3,k)*R(:,1)-s(1,k)*R(:,3) s(1,k)*R(:,2)-s(2,k)*R(:,1) I3]) )/ h(k)^2; 
end
initial = [0;0;0;t];
results = initial + (J' * J)\ (J') * (obs - f);
X_GN = results(1:3);
t_GN = results(4:6);
Xhat = [0 -X_GN(3) X_GN(2); X_GN(3) 0 -X_GN(1); -X_GN(2) X_GN(1) 0];
R_GN = R * expm(Xhat);   