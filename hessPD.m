function [hessP, hessD] = hessPD(A,B,C,beta,R,Rb)
% A, B, C, beta: M-by-1

singam = sqrt(Rb^2 - R^2*(cos(beta).^2))/Rb;
cosgam = R/Rb*cos(beta);
dgam = R/Rb*sin(beta)./singam;
ddgam = R/Rb*(singam.*cos(beta) - sin(beta).*cosgam.*dgam)./(singam.^2);
dC = R^2*sin(2*beta);
ddC = 2*R^2*cos(2*beta);
dsingam = cosgam.*ddgam - singam.*(dgam.^2);
tmp_BAC = B.^2 - 4*A.*C;

tmp_P0 = 2*Rb*A + B + 2*Rb*sqrt(A).*singam;
tmp_P1 = -4*A.*(tmp_BAC.*ddC + 4*A.*(dC.^2))./(tmp_BAC.^2);
tmp_P2 = -4*Rb*sign(tmp_P0).*sqrt(A).*(abs(tmp_P0).*dsingam - 2*Rb*sign(tmp_P0).*sqrt(A).*(cosgam.*dgam).^2)./(tmp_P0.^2);
hessP = 2*(Rb*(1-1./A).*dsingam + R*sin(beta) - B/4./(A.^(3/2)).*(tmp_P1 + tmp_P2));

tmp_D0 = singam + sqrt(C)/Rb + B./sqrt(C)/2;
tmp_lnD = log(tmp_BAC) - log(4*C) - 2*log(abs(tmp_D0));
tmp_dlnD0 = cosgam.*dgam + dC./sqrt(C)/2/Rb - B.*C.^(-3/2).*dC/4;
tmp_dlnD = -4*A.*dC./tmp_BAC - dC./C - 2*sign(tmp_D0).*tmp_dlnD0./abs(tmp_D0);
tmp_D1 = (C.^(3/2).*(dC.*sin(beta) + (2*C+ddC).*cos(beta)) - (2*C.*sin(beta) + cos(beta).*dC)*3/2.*sqrt(C).*dC)./(C.^3);
tmp_D2 = R/4./(C.^(3/2)).*(2*C.*sin(beta) + cos(beta).*dC);
tmp_D3 = -4*A.*(tmp_BAC.*ddC + 4*A.*(dC.^2))./(tmp_BAC.^2);
tmp_D4 = -(C.*ddC - dC.^2)./(C.^2);
tmp_D5 = -2*(tmp_D0.*(dsingam + (ddC./sqrt(C) - (dC.^2)./(C.^(3/2))/2)/2/Rb - B/4.*(ddC.*(C.^(-3/2)) - 3/2*(C.^(-5/2)).*(dC.^2))) - tmp_dlnD0.^2)./(tmp_D0.^2);
hessD = 2*R*(ddgam + R/4*tmp_D1.*tmp_lnD + 2*tmp_D2.*tmp_dlnD - R*cos(beta)./sqrt(C)/2.*(tmp_D3+tmp_D4+tmp_D5));

hessP = diag(hessP); hessD = diag(hessD);

end

