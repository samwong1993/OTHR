%% Zeroth-Order Algorithm for beta subprobem
%created by Huang Sen
%Email: huangsen1993@gmail.com
function beta = solve_eq(F,R,Rb,Rm,Ym,beta,XYZ,x,k)
ss = 0.001;
for i = 1:1000
[A B C] = ABC(F,R,Rb,Rm,Ym,beta);
[P D] = PD(A,B,C,beta,R,Rb);
dif_old = norm(XYZ(k,:)-x,2)^2 - 4*R^2*(1 - cos(D(k)/R))/2;
beta(k) = beta(k)-ss; 
[A B C] = ABC(F,R,Rb,Rm,Ym,beta);
[P D] = PD(A,B,C,beta,R,Rb);
dif = norm(XYZ(k,:)-x,2)^2 - 4*R^2*(1 - cos(D(k)/R))/2;
if sign(dif)~=sign(dif_old)
    ss = - ss/10;
else if abs(dif)>abs(dif_old)
    ss = -ss;
    end
end
if abs(dif)<1e-7
    break
end

end