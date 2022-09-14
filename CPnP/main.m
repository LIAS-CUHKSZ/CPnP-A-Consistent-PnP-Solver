% main codes: give an example to run the CPnP solver
% Copyright <2022>  <Guangyang Zeng, Shiyu Chen, Biqiang Mu, Guodong Shi, Junfeng Wu>
% Guangyang Zeng, SLAMLab-CUHKSZ, September 2022
% zengguangyang@cuhk.edu.cn, https://github.com/SLAMLab-CUHKSZ 
% paper link: https://arxiv.org/abs/2209.05824
clc
clear
rmse_t_1=zeros(1,7);
rmse_R_1=zeros(1,7);
rmse_t_GN=zeros(1,7);
rmse_R_GN=zeros(1,7);
CRB_R=zeros(1,7);
CRB_t=zeros(1,7);
monte_num=1;
counter=0;

for i=1:7
%% generate input data
N=round(10*10^((i-1)/(2)));
sigma=0.01;   
[A,point,Rt,focal]=generate_input_data(N,sigma);  % noise-free case
for j=1:monte_num
R=Rt(1:3,1:3);
t=Rt(1:3,4);
s=zeros(3,N);
Psens_2D=zeros(2,N);
for k=1:N
    s(:,k)=point(k).Xworld;
    Psens_2D(:,k)=(point(k).Ximg(1:2)+randn(2,1)*sigma*focal);  % add projection noise
end

%% estimate the pose and calculate rmse
[R_1,t_1,R_GN,t_GN]=CPnP(s,Psens_2D,focal,focal,A(1,3),A(2,3));

rmse_R_1(i)=rmse_R_1(i)+norm(R_1-R,"fro")^2;
rmse_t_1(i)=rmse_t_1(i)+norm(t-t_1)^2;
rmse_R_GN(i)=rmse_R_GN(i)+norm(R_GN-R,"fro")^2;
rmse_t_GN(i)=rmse_t_GN(i)+norm(t-t_GN)^2;
end
rmse_R_1(i)=sqrt(rmse_R_1(i)/monte_num);
rmse_t_1(i)=sqrt(rmse_t_1(i)/monte_num);
rmse_R_GN(i)=sqrt(rmse_R_GN(i)/monte_num);
rmse_t_GN(i)=sqrt(rmse_t_GN(i)/monte_num);

%% calculate the constrained CRB
F=zeros(12,12);
e1 = [1;0;0];
e2 = [0;1;0];
e3 = [0;0;1];
for k=1:N
    L = kron([s(:,k)' 1],eye(3));
    g1=e1'*L*[vec(R);t];
    g2=e2'*L*[vec(R);t];
    g3=e3'*L*[vec(R);t];
    F=F+((g3*L'*e1-g1*L'*e3)*(g3*e1'*L-g1*e3'*L)+(g3*L'*e2-g2*L'*e3)*(g3*e2'*L-g2*e3'*L))/(g3^4);
end
F=F/sigma^2;
O=zeros(3,1);
UU=([-R(:,3) O R(:,2) O O O;O -R(:,3) -R(:,1) O O O;R(:,1) R(:,2) O O O O;O O O e1*sqrt(2) e2*sqrt(2) e3*sqrt(2)])/sqrt(2);
bar_F=UU*((UU'*F*UU)\UU');
CRB_R(i)=sqrt(trace(bar_F(1:9,1:9)));
CRB_t(i)=sqrt(trace(bar_F(10:12,10:12)));
end
rmse_R_1
rmse_t_1
rmse_R_GN
rmse_t_GN
CRB_R
CRB_t