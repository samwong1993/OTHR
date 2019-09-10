%% Generate the max/min direct distance of feasible region
%created by Huang Sen
%Email: huangsen1993@gmail.com
function [max_dis,min_dis,upper] = beta_bound(M,F,R,Rb,Rm,Ym)
[A B C] = ABC(F,R,Rb,Rm,Ym,0);
[P D] = PD(A,B,C,0,R,Rb);
d = 2* R*sin(D/2/R);
max_dis = d;
if F<1
    con_beta = pi/2;
else
    A = 1-1/F^2+(Rb/F/Ym)^2;
    B = - 2*Rm*Rb^2/(F^2*Ym^2);
    B = B*ones(1,M);
    A = A*ones(1,M);
    con_beta = acos(sqrt(((Rb*Rm/F/Ym)^2.*ones(1,M) - B.^2/4./A)./R^2));
end
upper = UpBound_beta(F,R,Rb,Rm,Ym,con_beta(1));  
[A B C] = ABC(F,R,Rb,Rm,Ym,upper);
[P D] = PD(A,B,C,upper,R,Rb);
d = 2* R*sin(D/2/R);
min_dis = d;
end