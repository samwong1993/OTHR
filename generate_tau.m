%% Generate tau
%created by Huang Sen
%Email: huangsen1993@gmail.com
function [tau] = generate_tau(M,F,R,Rb,Rm,Ym,emitter,XYZ)
if F<1
    con_beta = pi/2;
else
    A = 1-1/F^2+(Rb/F/Ym)^2;
    B = - 2*Rm*Rb^2/(F^2*Ym^2);
    B = B*ones(1,M);
    A = A*ones(1,M);
    con_beta = acos(sqrt(((Rb*Rm/F/Ym)^2.*ones(1,M) - B.^2/4./A)./R^2));
end
beta = zeros(1,M);
[A B C] = ABC(F,R,Rb,Rm,Ym,beta);
[P D] = PD(A,B,C,beta,R,Rb);
x = emitter';
for i =1:7
    for k = 1:M
        beta = solve_eq(F,R,Rb,Rm,Ym,beta,XYZ,x,k);
    end
end
beta0 = beta;
[A B C] = ABC(F,R,Rb,Rm,Ym,beta0);
[P D] = PD(A,B,C,beta0,R,Rb);
%Generate tau
k = 1;
for i = 1:M
    for j = i+1:M
        tau(k) = P(i) - P(j);
        k = k + 1;
    end
end
end