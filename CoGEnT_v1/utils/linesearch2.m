function [ alpha1 ] = linesearch2( y,x )
%cvx_solver sdpt3
%cvx_begin quiet
%variable alpha1(1,1)
%minimize norm(y-x*alpha1)
%subject to 
%alpha1>=0
%alpha1<=1
%cvx_end
alpha1=x\y;
if alpha1<=0
    alpha1=0;
elseif alpha1>=1
    alpha1=1;
else
end

    

end

