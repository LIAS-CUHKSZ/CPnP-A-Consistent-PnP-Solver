% Copyright <2022>  <Guangyang Zeng, Shiyu Chen, Biqiang Mu, Guodong Shi, Junfeng Wu>
% Guangyang Zeng, SLAMLab-CUHKSZ, September 2022
% zengguangyang@cuhk.edu.cn, https://github.com/SLAMLab-CUHKSZ 
% paper link: https://arxiv.org/abs/2209.05824

function ratation_matrix=cal_rotation(a,b,r)
ratation_matrix=[cosd(a) -sind(a) 0;sind(a) cosd(a) 0;0 0 1]*[cosd(b) 0 sind(b);0 1 0;-sind(b) 0 cosd(b)]*[1 0 0;0 cosd(r) -sind(r);0 sind(r) cosd(r)];