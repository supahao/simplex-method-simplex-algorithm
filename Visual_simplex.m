function [x,fval,existFlag,BFS,A,b]=Visual_simplex(newf,NewAeq,Newbeq,B)
max_time_used = 3600;
tic;
%% 初始化解
[nrow,ncol]=size(NewAeq);
x=zeros(ncol,1);
fval=-1;
BFS = -2;
A = -1;
b = -1;
f=@(x) newf'*x;
existFlag=0;
it=0;
%% 单纯形法第一步，把当前表格变成规范型
% 如果不是规范型，then 先变成规范型
% 显式解直接套
% A=[NewAeq Newbeq];
%  for i = 1:nrow
%     if A(i,ncol+1)<-1e-7
%         A(i,:)=-A(i,:);
%     end
%  end
%  NewAeq=A(:,1:ncol);
%  Newbeq=A(:,ncol+1);
[newc,newAeq,newbeq,topright]=getCanonical(B,newf,NewAeq,Newbeq);

%% 得到规范型之后，判断当前基是否为最优解，不是的话做出基入基操作，转轴运算
c = [newc;topright];
A = [newAeq newbeq];
 while ~all(c(1:length(c)-1)>=-1e-7) 
  %如果做单纯形表法超时了，就返回超时标志
 disp("---------------------------------------------------------------------------------------------------------------------------");
 fprintf("迭代次数:%d\n当前单纯形表:\n",it);
 disp([c';A]);
 if toc>max_time_used
     existFlag=5;
     return;
 end
 d=find(c<-1e-7);%d(1)-------第一个负数元素列坐标
 e=find(A(:,d(1)) > 1e-7);% e包含d（1）列中正元素的行坐标，如果e没找到正元素，则说明该问题无界
 if isempty(e)
%      disp('该问题无界');
     existFlag=1;
     x=zeros(ncol,1);fval=-999;
     BFS=-1;b=-999;
     return
 end
 fprintf('x%d入基\n',d(1));
 g=[];
 for ii=1:length(e)
%  e是可以求比值的行index集合；遍历e中元素e(ii),此步用于寻找最大比率
 g=[g A(e(ii),ncol+1) / A(e(ii),d(1)) ];
 end
%  g中有d（1）对应的所有的比率，此处用最小比值法
 h=find(g==min(g));%选择离基变量,e(h)就是行坐标，离基变量的行index,也就是B中==e(h)的那一个要被替换成入基变量d(1)
if length(h)>=2
% 如果长度大于等于2，说明找到了相同的比值，用布兰德法则，选取下标最小的那个，也就是对应B中数字最小的那个
% 要比较e(h)对应的B中的值哪个小，选对应的最小的e(h)赋值给h
temp=max(B);
count=0;
for i = 1:length(h)
if B(e(h(i)))<temp
    temp = B(e(h(i)));
    count=i;
end
end
h = h(count);

end
 p=A(e(h),d(1));
fprintf("x%d离基\n",B(e(h)));
 B(e(h))=d(1);
 
%  对出基变量所在行操作，归一化
 A(e(h),:)=A(e(h),:)/p;
%  j是c应该加上的归一化后的出基行的倍数
 j=-c(d(1))/A(e(h),d(1));
 
%  对c进行操作
temprow = A(e(h),:);
 c = c' + j*temprow;
 c = c';
 
%  因为入基所在的列要符合单位阵的样式，所以其他行对应的列都要化为0
for ii=[1:e(h)-1,e(h)+1:nrow] %这个for就是把入基变量以外的行都归为0
%      j是要化为0的行对应元素的倍数
     j=-A(ii,d(1))/A(e(h),d(1));  
     A(ii,:) = j*A(e(h),:)+A(ii,:);
%    如果归零之后b中对应项变为负数，则整行*-1保证b全部>=0
     if A(ii,ncol+1)<0
        A(ii,:)=A(ii,:)*-1;
     end
 end%%%%%%%%%%%%截止，对A的操作完成
 it=it+1; %这个可以用来查看迭代的次数
 end
 disp('##########################################################################################################');
 disp('迭代完成时的单纯形表：');
 disp([c';A])
 
 
 %% 完成了寻找最优解的操作之后，只要用B将解还原出来即可
 BFS=B;
 [~,n]=size(A);
 b=A(:,n);
%  遍历B中元素，将x中对应的值赋为该行最后一个
 for i = 1:length(B)
 x(B(i)) = A(i,n);
 end
%  解已还原出来
fval=f(x);
%% 判断无穷解：如果有基对应的检验数=0就说明无穷解
rr = ones(nrow,1);
for i = 1:length(B)
    rr(i) = c(i);
end
if any(rr==0)
    existFlag=3;
end
    
    
    
end