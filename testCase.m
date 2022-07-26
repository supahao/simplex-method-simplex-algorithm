%% test simplex.m
% 用一个辅助问题做测试
clc;clear;
newf=[0;0;0;0;1;1];
NewAeq=[3 1 1 0 0 0 
        1 1 0 0 1 0 
        3 2 0 -1 0 1 ];
Newbeq=[27;12;30];    
B=[3,5,6];
[x,fval,existFlag,BFS,A,b]=simplex(newf,NewAeq,Newbeq,B)



%% findFirstBS test 1
clc;clear;
f=[1;1;1;0];
Aeq=[1 2 3 0
     -1 2 6 0
     0 4 9 0
     0 0 3 1];
beq=[3;2;5;1];
% A矩阵不是行满秩，要去掉多余行，变成行满秩
% How?
% 书上说，在返回的矩阵处理掉冗余的那一行然后删去那个基就行了
[optval,B,NewAeq,Newbeq] = findFirstBS(f,Aeq,beq)

%% test getCanonical
B=[3;5;6];
c=[0;0;0;0;1;1];
Aeq=[3 1 1 0 0 0
     1 1 0 0 1 0
     3 2 0 -1 0 1];
beq=[27;12;30];
[newc,newAeq,newbeq,topright]=getCanonical(B,c,Aeq,beq)
% 测试通过
%% getStandard Test1
% 测试第一种输入形式和求解情况
f=[1;-1];
A=[-1 -1 ;2 3];
b=[0 ; 8 ];
[f,A,b]=getStandard(f,A,b);
% [x,fval]=linprog(f,A,b)
% [x,fval,existFlag] = MyLPSolver(f,A,b);

%% getStandard Test2
% 测试第二种输入形式和求解情况
f=[3;-1;4];
A=[2 1 -4;2 6 1];
b=[8;0];
Aeq=[1 -1 -1; 0 1 2; 8 0 -1];
beq=[0 3 -1]';
[f,A,b]=getStandard(f,A,b,Aeq,beq)
% [x,fval]=linprog(f,A,b,Aeq,beq)
% [x,fval,existF] = MyLPSolver(f,A,b,Aeq,beq);

%% getStandard Test3
% 测试第三种输入形式和求解情况
f=[3;-1;4];
A=[2 1 -4;2 6 1];
b=[8;0];
Aeq=[1 -1 -1; 0 1 2; 8 0 -1];
beq=[0 3 -1]';
lb=[-1;-3;-2];
ub=[6 ; 6;6];
[f,A,b]=getStandard(f,A,b,Aeq,beq,lb,ub)
% [x,fval]=linprog(f,A,b,Aeq,beq,lb,ub)
% [x,fval,existFlag] = MyLPSolver(f,A,b,Aeq,beq,lb,ub);
% x = x + [lb;zeros(length(x)-length(lb),1)]%做这个变换，前面m个就是原解
% fval = fval + f'*lb %做这个加法就能得到原最优值

%% Generate_1 Test1 
%测试随机生成可行的线性规划问题函数
[f,A,b]=generate_1(300,170,0.99);
start = cputime;
[x,fval] = MyLPSolver(f,A,b);
time_used=cputime-start;
fprintf('用时为:%fs',time_used);

%% 作业中的手算问题1 测试
f=[-3;1;-2];
A=[1 3 1;-2 1 -1];
b=[5;-2];
Aeq=[4 3 -2];
beq=[5];
lb=[0;0;0];
ub=[9999;9999;9999];
% [~,linprog_f]=linprog(f,A,b,Aeq,beq,lb,ub)
[x,fval]=MyLPSolver(f,A,b,Aeq,beq,lb,ub)
%% 作业中的手算问题2 测试 本问题为无解问题
f=[1;-1;1];
A=[2 -1 -2;2 -3 -1;-1 1 1];
b=[4;-5;-1];
lb=[0;0;0];
ub=[9999;9999;9999];
[x,fval,existFlag]=MyLPSolver(f,A,b,[],[],lb,ub);

%% 输入维度不正确 测试1 
%形式1的测试
f=[1;-1];
A=[2 -1 -2;2 -3 -1;-1 1 1];
b=[4;-5;-1];
[x,fval,existFlag]=MyLPSolver(f,A,b);

%% 输入维度不正确 测试2
%形式2的测试
f=[-3;1;-2];
A=[1 3 1;3 3 3];
b=[5;-2];
Aeq=[4  -2];
beq=[];
[x,fval,existFlag]=MyLPSolver(f,A,b,Aeq,beq,[],[]);

%% 不满秩情况测试
%书本上的不满秩的简单情况用来做测试
f=[1;1;1;0];
Aeq=[1 2 3 0 
     -1 2 6 0
      0 4 9 0
      0 0 3 1];
beq=[3;2;5;1];
lb=[0;0;0;0];
ub=[99;99;99;99];
% [x,fv]=linprog(f,[],[],Aeq,beq,lb,ub)
[x,fval]=MyLPSolver(f,[],[],Aeq,beq,lb,ub)


%% 生成合理的不同规模的线性规划问题 测试2
[f,A,b]=generate_1(5,5,0.99,1);%最后还有一个flag：默认是1收敛，-3是无界，-2是无解
start = cputime;
[x,fval] = MyLPSolver(f,A,b);
time_used=cputime-start;
fprintf('用时为:%fs',time_used);

%% 可视化展示一个m=5,n=5的线性规划 
%要把MyLPSolver和findFirstBS中的simplex调用改为Visual_simplex!!!!!
f =[ 1
     4
     5
     1
     5];
A = [    0   -0.5342    0.2426         0         0
         0   -0.1006         0   -1.6250         0
   -1.5144         0    1.0262   -0.7581    2.0783
         0         0   -2.2220    0.4488         0
         0    0.0006   -0.7562    0.4043   -0.7939];
 b=[ 4
     2
    -9
    -9
    -7];
[x,fval]=MyLPSolver(f,A,b);

%% MPS文件输入求解 1
problem = mpsread("3kb-bk4x3.mps","ReturnNames",true);
problem.solver='linprog';
problem.options = optimoptions(problem.solver);
start = cputime;
[x,fval] = MyLPSolver(problem);
time_used=cputime-start;
x
fval
fprintf('用时为:%fs',time_used);
%% MPS文件输入求解 2
problem = mpsread("56kb-ran17x17.mps");
problem.solver='linprog';
problem.options = optimoptions(problem.solver);
start = cputime;
[~,fval] = MyLPSolver(problem);
time_used=cputime-start;
fprintf('用时为:%fs',time_used);
%% MPS文件输入求解 3
problem = mpsread("974kb-n3709.mps");
problem.solver='linprog';
problem.options = optimoptions(problem.solver);
start = cputime;
[~,fval] = MyLPSolver(problem);
time_used=cputime-start;
fval
fprintf('用时为:%fs',time_used);

%% m固定，n从小到大时，求解时间的变化
time_array=[];
n_array=[140:170];
for i=1:length(n_array)
[f,A,b]=generate_1(300,n_array(i),0.99,1);%最后还有一个flag：默认是1收敛，-3是无界，-2是无解
start = cputime;
[x,fval] = MyLPSolver(f,A,b);
time_used=cputime-start;
fprintf('用时为:%fs',time_used);
time_array(i) = time_used;
end
clf;
  title('m=300时，n从小到大求解时间的变化')
  xlabel('n的取值') 
  ylabel('求解时间(s)')
  hold on
plot(n_array,time_array,'-x');

%% A是否为稀疏矩阵时 求解时间的变化
%同一m，n下，多次的不同稀疏度的耗时情况
time_array=zeros(2,10);
n_array=[140:170];
m=400;
n=200;
count =10;
while count ~=0

[f,A,b]=generate_1(m,n,0.999,1);%最后还有一个flag：默认是1收敛，-3是无界，-2是无解
start = cputime;
[x,fval] = MyLPSolver(f,A,b);
time_used=cputime-start;
fprintf('用时为:%fs',time_used);
time_array(1,count) = time_used;

[f,A,b]=generate_1(m,n,0.05,1);%最后还有一个flag：默认是1收敛，-3是无界，-2是无解
start = cputime;
[x,fval] = MyLPSolver(f,A,b);
time_used=cputime-start;
fprintf('用时为:%fs',time_used);
time_array(2,count) = time_used;
count = count-1;

end
clf;
  title('m=400,n=200，稠密与稀疏的对比')
  xlabel('左为稠密，右为稀疏') 
  ylabel('求解时间(s)')
  hold on
bar(time_array);




%% MPS文件输入求解 4
problem = mpsread("20kb-bal.mps");
problem.solver='linprog';
problem.options = optimoptions(problem.solver);
start = cputime;
[~,fval] = MyLPSolver(problem);
time_used=cputime-start;
fval
fprintf('用时为:%fs',time_used);




