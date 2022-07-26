function [newf,newAeq,newbeq]=getStandard(f,A,b,Aeq,beq,lb,ub)
[m,n]=size(A);

if nargin==3
    newf=[f;-f];
    newAeq=[A -A];
    newf = [newf ; zeros(m,1)];
    newAeq=[newAeq eye(m)];
    newbeq=b;
    
elseif nargin==5
    [p,q]=size(Aeq);
    newf=[f;-f;zeros(m,1)];
    % newA=[A -A eye(m)];
    % newAeq=[Aeq -Aeq zeros(p,m)];
    % newAeq = [newA;newAeq]
    newAeq =[A -A eye(m);Aeq -Aeq zeros(p,m)];
    newbeq=[b;beq]
    
elseif nargin==7
    [p,q]=size(Aeq);
    if n==0
        n = q;
    end
    [p,q]=size(Aeq);
    if ~isempty(b)
        b = b-A*lb;
    end
    if ~isempty(beq)
        beq = beq-Aeq*lb;
    end
    newf=[f;zeros(m+length(f),1)];
    ub=ub-lb;
    A = [A eye(m)];
    Aeq = [Aeq zeros(p,m)];
    newAeq = [A;Aeq];
    newAeq = [newAeq zeros(m+p,n)];
    newAeq = [newAeq;eye(n) zeros(n,m) eye(n)];
    newbeq = [b;beq;ub];
else
    disp('输入不符合要求，请按规定的形式进行输入')
end

end