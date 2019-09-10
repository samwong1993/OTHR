%% Calculate the gradient of objective function
%created by Huang Sen
%Email: huangsen1993@gmail.com
function graF = GraF(A,B,C,beta,R,Rb,inv,tau,delta,G)
[graP graD]=graPD(A,B,C,beta,R,Rb);
[P D] = PD(A,B,C,beta,R,Rb);
gra1 = (2*G'*inv*(G*P'-tau') + delta*ones(length(P),1))'.*graP; %- 2*R^2*sin(D/R)/R.*graD.*lambda'; %- rho*penalty'*2*R.*sin(D/R)/R.*graD;
gra2 = 0;%(lambda)'.*(- 2*R^2*sin(D/R)/R.*graD);
graF = gra1 + gra2;
end