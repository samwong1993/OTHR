%% Calculate the upper bound of beta/The minimum for D
%created by Huang Sen
%Email: huangsen1993@gmail.com
function beta = CD(F,R,Rb,Rm,Ym,upper,G,tau,inv_Omega,beta,k)
beta4 = beta;
beta1 = beta;
beta2 = beta;
beta3 = beta;
beta4(k) = upper;
beta1(k) = 0;
r = 0.618;
for i = 1:5000
beta2(k) = r*beta1(k) + (1- r)*beta4(k);
beta3(k) = r*beta4(k) + (1- r)*beta1(k);
[A B C] = ABC(F,R,Rb,Rm,Ym,beta1);
[P D] = PD(A,B,C,beta1,R,Rb);
D1 = (G*P'-tau')'*inv_Omega*(G*P'-tau');
[A B C] = ABC(F,R,Rb,Rm,Ym,beta2);
[P D] = PD(A,B,C,beta2,R,Rb);
D2 = (G*P'-tau')'*inv_Omega*(G*P'-tau');
[A B C] = ABC(F,R,Rb,Rm,Ym,beta3);
[P D] = PD(A,B,C,beta3,R,Rb);
D3 = (G*P'-tau')'*inv_Omega*(G*P'-tau');
[A B C] = ABC(F,R,Rb,Rm,Ym,beta4);
[P D] = PD(A,B,C,beta4,R,Rb);
D4 = (G*P'-tau')'*inv_Omega*(G*P'-tau');
if D2<D3
    beta4(k) = beta3(k);
else
    beta1(k) = beta2(k);
end
if abs(D1-D2)<1e-7
    break
end
end
beta(k) = (beta1(k)+beta2(k))/2;
end