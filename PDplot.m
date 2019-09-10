%% Plot the function value of P/D
%created by Huang Sen
%Email: huangsen1993@gmail.com
clear
R = 6371.2;
M = 5;
N = M*(M-1)/2;
Omega = 0.5*(ones(N,N) + eye(N));
inv_Omega = Omega^-1;
delta = 0;
Rm = 6650;
Ym = 100;
Rb = Rm - Ym;
fc = 10;
f = 1;
F = f/fc;
if F<1
    con_beta = pi/2;
else
    A = 1-1/F^2+(Rb/F/Ym)^2;
    B = - 2*Rm*Rb^2/(F^2*Ym^2);
    B = B*ones(1,M);
    A = A*ones(1,M);
    con_beta = acos(sqrt(((Rb*Rm/F/Ym)^2.*ones(1,M) - B.^2/4./A)./R^2));
end




beta = 0:0.001:pi;
[graP graD] = graPD(A,B,C,beta,R,Rb);
[A B C] = ABC(F,R,Rb,Rm,Ym,beta);
[P D] = PD(A,B,C,beta,R,Rb);
r = -B/2./A;
subplot(1,2,1)
plot(beta,D, 'k-', 'linewidth', 1.1)
xlabel('beta')
ylabel('D')
subplot(1,2,2)
plot(beta,P,'k-', 'linewidth', 1.1)
xlabel('beta')
ylabel('P')