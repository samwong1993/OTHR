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
f = 15;
F = f/fc;
%perturbation on ionosphere
P_Rm = Rm;
P_Ym = Ym;
P_Rb = P_Rm - P_Ym;
P_fc = fc;
P_f = f;
P_F = P_f/P_fc;
[max_dis,min_dis,upper] = beta_bound(M,F,R,Rb,Rm,Ym);
[G] = generate_G(N,M);
%[emitter,XYZ,beta0] = generator(M,F,R,Rb,Rm,Ym,max_dis,min_dis);
beta0 = [0.114957231412252,0.449398124172348,0.277420425918117,0.0168095219080640,0.103488345084960];
XYZ = zeros(M,3);
%Hong Kong
[x0 y0 z0] = LGLTtoXYZ(114.16,22.28,R);
emitter = [x0 y0 z0]';
%Bei Jing
[x0 y0 z0] = LGLTtoXYZ(116.41,39.90,R);
XYZ(1,:) = [x0 y0 z0];
%Wu Han
[x0 y0 z0] = LGLTtoXYZ(114.31,30.59,R);
XYZ(2,:) = [x0 y0 z0];
%Shang Hai
[x0 y0 z0] = LGLTtoXYZ(121.47,31.23,R);
XYZ(3,:) = [x0 y0 z0];
%Tokyo
[x0 y0 z0] = LGLTtoXYZ(139.69,35.69,R);
XYZ(4,:) = [x0 y0 z0];
%Seoul
[x0 y0 z0] = LGLTtoXYZ(126.58,37.33,R);
XYZ(5,:) = [x0 y0 z0];
Omega = 0.5*(ones(N,N) + eye(N));
inv_Omega = Omega^-1;
tau = generate_tau(M,F,R,Rb,Rm,Ym,emitter,XYZ);
%Step 3
ss = 1e-9;
x = [];
x = Initialization(XYZ,max_dis,min_dis,M,inv_Omega,G,tau,F,R,Rb,Rm,Ym,upper,x);
beta = zeros(1,M);
[A B C] = ABC(F,R,Rb,Rm,Ym,beta);
[P, D] = PD(A,B,C,beta,R,Rb);
obj_old = (G*P'-tau')'*inv_Omega*(G*P'-tau');
[graP, graD] = graPD(A,B,C,beta,R,Rb);
dP_beta = (2*G'*inv_Omega*(G*P'-tau'))'.*graP;
beta = beta - ss*dP_beta;
beta(beta<0) = 0;
beta(beta>upper) = upper;
[A B C] = ABC(F,R,Rb,Rm,Ym,beta);
[P, D] = PD(A,B,C,beta,R,Rb);
obj = (G*P'-tau')'*inv_Omega*(G*P'-tau');
for i = 1:1000
    [A B C] = ABC(F,R,Rb,Rm,Ym,beta);
    [P, D] = PD(A,B,C,beta,R,Rb);
    r = 2*R*sin(D./2/R);
    for j = 1:M 
        x = (x - XYZ(j,:))/norm(x - XYZ(j,:))*r(j)+XYZ(j,:);
    end
    x = x/norm(x)*R;
end
for k = 1:M
    beta_old = solve_eq(F,R,Rb,Rm,Ym,beta,XYZ,x,k);
end
if norm(beta_old-beta)<8e-4
    fprintf("success\n")
end
obj
