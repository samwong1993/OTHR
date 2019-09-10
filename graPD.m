function [graP, graD] = graPD(A,B,C,beta,R,Rb)

dgam = R*sin(beta)./sqrt(Rb^2 - R^2*(cos(beta).^2));
singam = sqrt(Rb^2 - R^2*(cos(beta).^2))/Rb;
cosgam = R/Rb*cos(beta);
dC = R^2*sin(2*beta);

tmp_denom = 2*Rb*A + B + 2*Rb*sqrt(A).*singam;
graP = 2*(Rb*(1 - 1./A).*cosgam.*dgam - R*cos(beta) - B/4./A.^(3/2).*(-4*A.*dC./((B.^2 - 4.*A.*C)) - 4*sign(tmp_denom)*Rb.*sqrt(A).*cosgam.*dgam./abs(tmp_denom)));
% graP = 2*R*((1-1./A).*cos(beta).*dgam - cos(beta) + B./sqrt(A)*R*sin(2*beta)./(B.^2 - 4.*A.*C) + B./A.*cos(beta).*dgam./tmp_denom);

tmp2_lnD = singam + sqrt(C)/Rb +B./sqrt(C)/2;
tmp_lnD = log(B.^2 - 4*A.*C) - log(4*C) - 2*log(abs(tmp2_lnD));
tmp_dlnD = 2*sign(tmp2_lnD).*(cosgam.*dgam + dC./sqrt(C)/2/Rb - B.*C.^(-3/2).*dC/4)./abs(tmp2_lnD);
graD = 2*R*(dgam - 1 - R/2./C.*(-sqrt(C).*sin(beta) - cos(beta)./sqrt(C).*dC/2).*tmp_lnD - R*cos(beta)./sqrt(C)/2.*(-4*A.*dC./(B.^2 - 4*A.*C) - dC./C - tmp_dlnD));

end

