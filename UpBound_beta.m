%% Calculate the upper bound of beta/The minimum for D
%created by Huang Sen
%Email: huangsen1993@gmail.com
function beta = UpBound_beta(F,R,Rb,Rm,Ym,con)
beta4 = con-1e-5;
beta1 = 0;
r = 0.618;
for i = 1:5000
beta2 = r*beta1 + (1- r)*beta4;
beta3 = r*beta4 + (1- r)*beta1;
[A B C] = ABC(F,R,Rb,Rm,Ym,beta1);
[P D] = PD(A,B,C,beta1,R,Rb);
D1 = D;
[A B C] = ABC(F,R,Rb,Rm,Ym,beta2);
[P D] = PD(A,B,C,beta2,R,Rb);
D2 = D;
[A B C] = ABC(F,R,Rb,Rm,Ym,beta3);
[P D] = PD(A,B,C,beta3,R,Rb);
D3 = D;
[A B C] = ABC(F,R,Rb,Rm,Ym,beta4);
[P D] = PD(A,B,C,beta4,R,Rb);
D4 = D;
if D2<D3
    beta4 = beta3;
else
    beta1 = beta2;
end
if abs(D1-D2)<1e-7
    break
end
end
beta = (beta1+beta2)/2;
end
