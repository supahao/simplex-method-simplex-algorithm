function [f,A,b]=generate_1(m,n,K,flag)
%% 超参数初始化
if nargin ==3
    flag = 1
end
% n=200;%变量个数
%目标函数系数的上下限
f_lb=-5;
f_ub=5;

%矩阵大小
% row = 400;
col = n;
% K=0.1;  % sparseness #非零元素占比
A_lb = -50;
A_ub = 50;

%b的设置
b_lb=-100;
b_ub=100;
%% 生成随机线性规划
f = randi([f_lb,f_ub],n,1);
% A = randi([A_lb,A_ub],row,col);
A = sprandn(m,col,K)*A_ub;
A = full(A);
b = randi([b_lb,b_ub],m,1);

[x,fval,flag_F] = linprog(f,A,b)    ;
if flag_F~=flag
    [f,A,b]=generate_1(m,n,K,flag);
    

end