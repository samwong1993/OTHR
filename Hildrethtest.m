clear
R = 6371.2;
M = 3;
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
% %Tokyo
% [x0 y0 z0] = LGLTtoXYZ(139.69,35.69,R);
% XYZ(4,:) = [x0 y0 z0];
% %Seoul
% [x0 y0 z0] = LGLTtoXYZ(126.58,37.33,R);
% XYZ(5,:) = [x0 y0 z0];
x_opt = randn(3,1);
x_opt = x_opt/norm(x_opt)*R;
x = x_opt;
b = zeros(M,1);
for i = 1:M
    A(2*i-1,:) = XYZ(i,:);
    b(2*i-1) = R^2 - min_dis^2/2;
    A(2*i,:) = - XYZ(i,:);
    b(2*i) = max_dis^2/2 - R^2;
end
A*x<=b


cvx_begin
    variable x(3)
    minimize( norm(x-x_opt) )
    subject to
        R^2 - max_dis^2/2<= XYZ(1,:)*x<= R^2 - min_dis^2/2
        R^2 - max_dis^2/2<= XYZ(2,:)*x<= R^2 - min_dis^2/2
        R^2 - max_dis^2/2<= XYZ(3,:)*x<= R^2 - min_dis^2/2
cvx_end
dis_old = norm(x-x_opt);
A*x<=b

x = proj_poly(XYZ,max_dis,min_dis,M,R,x);
x = x/norm(x)*R;
x
dis = [];
for iter = 1:M
    dis = [dis norm(XYZ(iter,:) - x')];
end
all(dis<max_dis)&all(dis>min_dis)