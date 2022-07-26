function [optval,BFS,NewAeq,Newbeq,existFlag]=findFirstBS(f,Aeq,beq)
% % 已经构造好标准型，寻找初始可行基需要构造辅助问题

%:nrow:约束数量
%:ncol:变量数量
[nrow,ncol]=size(Aeq);

%% 辅助问题的目标函数向量长度为 原变量数+约束数量
newf=zeros(ncol+nrow,1);
for i = ncol+1:nrow+ncol
    newf(i,1)=1;
end



%% 新的Aeq需要在右侧加上单位阵
NewAeq = [Aeq,eye(nrow)];
% 添加人工变量，i.e. 构造好辅助问题之后就要求解辅助问题

%% 有了辅助问题，就要用单纯形法求解
B=find(newf~=0);
Newbeq=beq;

%% 用单纯形法求解辅助问题 
[~,fval,existFlag,BFS,A,b] = simplex(newf,NewAeq,Newbeq,B);
% [~,fval,existFlag,BFS,A,b] = Visual_simplex(newf,NewAeq,Newbeq,B);
if existFlag == 5 
    optval=-1;
    return
end
if BFS==-1
    disp('求解辅助问题时出现无界情况');
end
%% 单纯形法求解好了辅助问题之后，如果fval=0则有初始可行基；
% 否则fval>0说明无解。无论怎样，本函数是要返回 最优值、初始可行基、已经完成部分计算的NewAeq和Newbeq（optval,BFS,NewAeq,Newbeq）
% NewAeq要去掉A中构建的人工变量对应的全部的列
optval=fval;
NewAeq = A( : , 1: ncol ) ;
Newbeq=b;


end