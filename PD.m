%% Calculate the function value of P/D
%created by Huang Sen
%Email: huangsen1993@gmail.com
function [P D] = PD(A,B,C,beta,R,Rb)
Fai = acos( R/Rb*cos(beta));
D = 2*R*( (Fai - beta) - R.*cos(beta)/2./sqrt(C).*log((B.^2-4*A.*C)./(4.*C.*(sin(Fai)+1/Rb*sqrt(C)+1/2./sqrt(C).*B).^2)) );
P = 2*( Rb*sin(Fai) - R*sin(beta) + 1./A.*(- Rb*sin(Fai)- B/4./sqrt(A).*log( (B.^2-4.*A.*C)./(2*A*Rb+B+2*Rb*sqrt(A).*sin(Fai)).^2) ));
end