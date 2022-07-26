function [x,fval,existFlag]=MyLPSolver(f,A,b,Aeq,beq,lb,ub)

%初始化解
x=-1; %最优解
fval=-1; %最优值
existFlag = -1; %结果标志
flag = -1; %输入形式标志


%检测是否是MPS输入
problemInput = false;
if nargin == 1
    if isa(f,'struct')
        problemInput = true;
        [f,A,b,Aeq,beq,lb,ub] = separate_Struct(f);
%        mpsread 将 problem.Aineq 和 problem.Aeq 以稀疏矩阵形式返回，转为完整矩阵
        A = full(A);Aeq=full(Aeq);
    else % Single input and non-structure.
        error(message('optim:linprog:InputArg'));
    end
    flag = 3;
    [newf,newAeq,newbeq]=getStandard(f,A,b,Aeq,beq,lb,ub);
end

[m,n] = size(A);
if n==0 && ~isempty(Aeq)
    [p,q]=size(Aeq);
    n=q;
end
%先转化成标准型然后两步法求解
if nargin==3 %第一种输入形式
    flag = 1;
    if length(f)~=n
       existFlag = 4;
       disp("f长度与A矩阵列数不匹配");
       return;
    end
    if length(b)~=m
        existFlag = 4;
        disp("A矩阵行数与b长度不匹配");
        return
    end
    [newf,newAeq,newbeq]=getStandard(f,A,b);
end
if nargin==5 %第二种输入形式
    [p,q]=size(Aeq);
    flag = 2;
    if ~isempty(A)
    if length(f)~=n  
       existFlag = 4;
       disp("f长度与A矩阵列数不匹配");
       return;
    end
    end
    if length(b)~=m
        existFlag = 4;
        disp("A矩阵行数与b长度不匹配");
        return
    end
    if length(f)~=q
        existFlag=4;
        disp("f的长度和Aeq的列数不相等");
        return;
    end
    if length(beq)~=p
        existFlag=4;
        disp("beq的长度和Aeq的行数不相等");
        return;
    end
    
    [newf,newAeq,newbeq]=getStandard(f,A,b,Aeq,beq);

end
if nargin==7 %第三种输入形式
    if m==0
        
    elseif m~=length(b)
        disp('A的大小与b不匹配');
        return
    end
    flag = 3;
    [newf,newAeq,newbeq]=getStandard(f,A,b,Aeq,beq,lb,ub);
end

%得到标准型之后，要保证newbeq如果有Inf，转换为可以求解的值，用9999代替Inf
find_Inf = find(newbeq==Inf);
for i = 1:length(find_Inf)
   newbeq(find_Inf(i)) = 99; 
end


%得到标准型后，要保证newbeq中都非负，这是两步法的起点，必须处理一下
for i = 1:length(newbeq)
    if newbeq(i)<0
        newbeq(i) = -newbeq(i);
        newAeq(i,:) = -newAeq(i,:);
    end
end

[new_nrow,new_ncol] = size(newAeq);

rankA = rank(newAeq);

if rankA~=new_nrow
    disp('标准化后的Aeq不满秩');
    %不满秩要高斯消元找出冗余行并删去它
    %此处判断矛盾约束导致无解:拼起来之前的newAeq不满秩，拼好之后满秩则说明有矛盾约束，直接无解
    tempAeq=[newAeq newbeq];
    rankTempAeq=rank(tempAeq);
    gap = new_nrow-rankA;
    if rankA~=rankTempAeq
        existFlag = 2;
        error("约束中存在矛盾，无解");
    end
    %高斯消元有点难做，不如一行行消除然后检查rank，如果删去了某一行后newAeq矩阵的秩不减，说明这一行没有用
       while gap~=0
        temp_row=size(newAeq);
        for i = temp_row:-1:1
            %遍历每一行
            tempAeq = newAeq;
            tempAeq(i,:)=[];
            if rank(tempAeq)==rankA
                newAeq(i,:)=[];
                newbeq(i)=[];
                gap = gap-1;
                break;
            end
        end
        
    end
    disp('已将冗余行去除');
end
% 化为标准型后，要用两步法求解

% 两步法-第一步，寻找初始可行基，找不到说明问题无解
[optval,B,NewAeq,Newbeq,existFlag]=findFirstBS(newf,newAeq,newbeq);

if existFlag~=5
    
    
    if optval~=0
        % 问题无解
        existFlag=2;
    else
        
        % 两步法-第二步，在初始可行基的基础上，化为规范型，继续做单纯形法直到停止
        % 可能找到唯一最优解、无数最优解、可能问题无界、
        % 用newf、NewAeq、Newbeq、B继续做单纯形法
        % 需要的返回结果是x,fval,existFlag：最优解，最优值，状态标志
        [x,fval,existFlag] = simplex(newf,NewAeq,Newbeq,B);
%         [x,fval,existFlag] = Visual_simplex(newf,NewAeq,Newbeq,B);
        
    end
    
end
    
    switch existFlag
        case 0
            %只有单个解
            disp('找到最优解，有1个');
            [x,fval]=recover_Solu(flag,f,lb,fval,x,n);
        case 1
            %线性规划无界
            disp('线性规划无界');
            return;
        case 2
            %线性规划无解
            disp('线性规划无解');
            return;
        case 3
            %有无数个最优解
            disp('找到最优解，有无数个,给出其中一个');
            if flag~=3
            [x,fval]=recover_Solu(flag,f,fval,x,n);
            else
            [x,fval]=recover_Solu(flag,f,fval,x,n,lb);
            end
        case 4
            %输入形式错误
            disp('输入形式错误')
        case 5
            %求解超时
            disp('求解超时')
    end
    

end


function [x,fval]=recover_Solu(flag,f,fval,x,n,lb)

switch flag
    case 1
        %第一种情况，恢复原解，原解的长度为n
        temp = zeros(n,1);
        for i = 1:n
            temp(i) = x(i)-x(i+n);
        end
        x = temp;
    case 2
        %第二种情况，恢复原解，原解的长度为n
        temp = zeros(n,1);
        for i = 1:n
            temp(i) = x(i)-x(i+n);
        end
        x = temp;
        
    case 3
        %第三种输入，恢复原解和最优值,原解的长度为n
        x = x + [lb ; zeros( length(x)-length(lb) ,1)];
        x = x(1:n);
        fval = f'*x;
end
end