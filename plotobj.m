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
beta0 = [0.349025176895742,0.279243732945529,0.174632834143026,0.183395025597251,0.285385790728298];
%beta0 = [0.674175558279844,0.678354140357644,0.679512211766613,0.679496774102502,0.678155318189844];
XYZ = zeros(M,3);
[x0 y0 z0] = LGLTtoXYZ(116.24,39.55,R);
emitter = [x0 y0 z0]';
[x0 y0 z0] = LGLTtoXYZ(128.72,40.55,R);
XYZ(1,:) = [x0 y0 z0];
[x0 y0 z0] = LGLTtoXYZ(130.42,38.68,R);
XYZ(2,:) = [x0 y0 z0];
[x0 y0 z0] = LGLTtoXYZ(132.94,33.82,R);
XYZ(3,:) = [x0 y0 z0];
[x0 y0 z0] = LGLTtoXYZ(130.90,31.84,R);
XYZ(4,:) = [x0 y0 z0];
[x0 y0 z0] = LGLTtoXYZ(129.06,35.63,R);
XYZ(5,:) = [x0 y0 z0];
[max_dis,min_dis,upper] = beta_bound(M,F,R,Rb,Rm,Ym);
[G] = generate_G(N,M);
tau = generate_tau(M,F,R,Rb,Rm,Ym,emitter,XYZ);
num = 200;
[x0, y0, z0] = ellipsoid(0, 0, 0, R, R, R, num);
for i = 1:num+1
    for j = 1:num+1
        x = [x0(i,j) y0(i,j) z0(i,j)];
        dis = [];
        for iter = 1:M
            dis = [dis norm(XYZ(iter,:) - x)];
        end
        %dis = [norm(XYZ(1,:) - x) norm(XYZ(2,:) - x) norm(XYZ(3,:) - x) norm(XYZ(4,:) - x) norm(XYZ(5,:) - x)];
        if all(dis<max_dis)&all(dis>min_dis)
            beta = 0.1*ones(1,M);
            for iter = 1:20
                for k = 1:M
                    beta = solve_eq(F,R,Rb,Rm,Ym,beta,XYZ,x,k);
                end
            end
            if all(beta>0)
                [A B C] = ABC(F,R,Rb,Rm,Ym,beta);
                [P D] = PD(A,B,C,beta,R,Rb);
                obj(i,j) = (G*P'-tau')'*inv_Omega*(G*P'-tau');
            else
                obj(i,j) = 0;
            end
        else
            obj(i,j) = 0;
        end
    end
end
% obj=obj/max(max(obj));
obj(obj==1e7)=5e6
colormap('hot')
imagesc(obj)
surf(obj)
xlim([1 200])
ylim([1 200])
% shading interp
axis off

obj_min = 9999999;
for i = 1:num+1
    for j = 1:num+1
        x = [x0(i,j) y0(i,j) z0(i,j)];
        dis = norm(emitter' - x);
        if dis<obj_min
            obj_min = dis;
            m = i;
            n = j;
        end
    end
end