function [newc,newAeq,newbeq,topright]=getCanonical(B,c,Aeq,beq)
% 获得了标准型和可行基，需要将标准型转化为规范型，然后传回f,c,Aeq,b回去
% 已有标准型和可行基要转为该基下的标准型只要做求逆运算即可，套公式
A_B=[];
c_B=[];
for i=1:length(B)
    A_B=[A_B Aeq(:,B(i))];
    c_B=[c_B;c(B(i))];    
end
inv_A_B=inv(A_B);
newAeq=inv_A_B*Aeq;
newbeq=inv_A_B*beq;
newc=c'-c_B'*inv_A_B*Aeq;
newc=newc';
topright=-c_B'*inv_A_B*beq;
end