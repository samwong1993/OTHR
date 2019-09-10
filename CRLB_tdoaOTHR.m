function CRLB = CRLB_tdoaOTHR(F, Rb, Ym, Rm, R, beta0, XYZ0, emitter, sigma, M,index)
% remark: sigma is a scalar
% betaTrue (M-by-1) is the true flying angle
% SensorPositions (M-by-3)
% CRLB_tdoaOTHR_new(F, Rb, Ym, Rm, R, beta0, XYZ, emitter, sigma, M)
for i = 1:M
    SensorPositions(i,:) = XYZ0(index(i),:);
    betaTrue(i) = beta0(index(i));
end
[A, B, C] = ABC(F,R,Rb,Rm,Ym,betaTrue);
[graP, graD] = graPD(A,B,C,betaTrue,R,Rb);

N = M*(M-1)/2;
dP = zeros(N,M); tmp_cnt = 0;
[G] = generate_G(N,M);
dP = G*diag(graP);
% for i = 1:M-1
%     for j = i+1:M
%         dP(tmp_cnt+(j-i),i) = graP(i);
%         dP(tmp_cnt+(j-i),j) = -graP(j);
%     end
%     tmp_cnt = tmp_cnt + M-i;
% end

U = 0.5*(ones(N,N) + eye(N));
invR = sigma^(-2)*U^-1;

[~, D] = PD(A,B,C,betaTrue,R,Rb);
S = SensorPositions; x = emitter;
dBeta = zeros(M, 2);
for i = 1:M
    dBeta(i,1) = -(S(i,1) - S(i,3)*x(1)/x(3))/R/sin(D(i)/R)/graD(i);
    dBeta(i,2) = -(S(i,2) - S(i,3)*x(2)/x(3))/R/sin(D(i)/R)/graD(i);
end

J = dBeta'*dP'*invR*dP*dBeta; %FIM
CRLB = sqrt(sum(diag(inv(J)))); %MSE
end

