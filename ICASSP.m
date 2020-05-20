clear
fid=fopen('randM5WTtau0_0001.txt','a+');
R = 6371.2;
%plt = 1 plot the earth,emitter,sensors and the generated sequence
plt = 0;
if plt == 1
    figure('color','k')
    %Add legend
    %emitter
    point1 = scatter3(0,0,0,50,'filled','r');
    hold on
    %sensors
    point2 = scatter3(0,0,0,'filled','c');
    %initial point
    point3 = scatter3(0,0,0,50,'filled','b');
    %recover location
    point4 = scatter3(0,0,0,50,'k*');
    %generated sequence
    point5 = scatter3(0,0,0,5,'m');
    %plot earth
    npanels = 72;
    alpha   = 1;
    image_file = 'earth.jpg';
    [x0, y0, z0] = ellipsoid(0, 0, 0, R, R, R);
    globe = surf(x0, y0, -z0, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
    cdata = imread(image_file);
    set(gca, 'NextPlot','add', 'Visible','off');
    set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');
    axis off;
    axis equal;
    axis auto;
end
M = 5;
N = M*(M-1)/2;
Omega = covariance(1,M);
inv_Omega = inv(Omega(1:M-1,1:M-1));
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
delta = 0.0001;
[max_dis,min_dis,upper] = beta_bound(M,F,R,Rb,Rm,Ym);
[G] = generate_G(N,M);
%[emitter,XYZ,beta0] = generator(M,F,R,Rb,Rm,Ym,max_dis,min_dis);
% beta0 = [0.277420425918117,0.0168095219080640,0.103488345084960];
% XYZ = zeros(M,3);
% %Hong Kong
% [x0 y0 z0] = LGLTtoXYZ(114.16,22.28,R);
% emitter = [x0 y0 z0]';
% %Shang Hai
% [x0 y0 z0] = LGLTtoXYZ(121.47,31.23,R);
% XYZ(1,:) = [x0 y0 z0];
% %Tokyo
% [x0 y0 z0] = LGLTtoXYZ(139.69,35.69,R);
% XYZ(2,:) = [x0 y0 z0];
% %Seoul
% [x0 y0 z0] = LGLTtoXYZ(126.58,37.33,R);
% XYZ(3,:) = [x0 y0 z0];
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

% choose = [2 4 5]
% beta0 = beta0(choose);
% XYZ = XYZ(choose,:);
sigma = [0:100:1000];
for index = 22:100%1:length(alpha)
    for noise_level = 1:length(sigma)
        sigma_t = sigma(noise_level)* 10^-9 * 3 * 10^5 ;
        [G] = generate_G(N,M);
        noise_t0 = randn(M,1);
        noise_t = (sigma_t*G*noise_t0)';
        tau = generate_tau(M,F,R,Rb,Rm,Ym,emitter,XYZ) + noise_t;
        G = G(1:M-1,:);
        tau = tau(1:M-1);
        %Step 1 Coordernate Descent
        beta = zeros(1,M);
        [A B C] = ABC(F,R,Rb,Rm,Ym,beta);
        [P D] = PD(A,B,C,beta,R,Rb);
        obj_old = (G*P'-tau')'*inv_Omega*(G*P'-tau') + delta*ones(M,1)'*P';
        for i =1:20
            for k = 1:M
                beta = CD(F,R,Rb,Rm,Ym,upper,G,tau,inv_Omega,beta,k);
                [A B C] = ABC(F,R,Rb,Rm,Ym,beta);
                [P D] = PD(A,B,C,beta,R,Rb);
            end
        end
        %Step 2 Inialization
        [A B C] = ABC(F,R,Rb,Rm,Ym,beta);
        [P D] = PD(A,B,C,beta,R,Rb);
        Xi = [];
        g = [];
        S = [];
        for i = 1:M
            Xi = [Xi XYZ(i,:)'];
            S = [S; norm(XYZ(i,:))^2+R^2];
            g = [g;(2*R*sin(D(i)/2/R))^2];
        end
        Xi = 2*Xi';
        A_xi = (Xi'*Xi)^-1*Xi';
        x0 = A_xi*(S - g);
        norm(x0 - emitter)
        x = x0';
        for i =1:1
            for k = 1:M
                beta = solve_eq(F,R,Rb,Rm,Ym,beta,XYZ,x,k);
            end
        end
        % beta(beta<0) = 0;
        % beta(beta>upper) = upper;
        x0 = x;
        beta = zeros(1,M);
        for i =1:1
            for k = 1:M
                beta = solve_eq(F,R,Rb,Rm,Ym,beta,XYZ,x,k);
            end
        end
        %Step 3
        thres = 10;
        alpha = 1.05;
        ss = 1e-4;
        center = R*sum(XYZ)/M/(norm(sum(XYZ)/M));
        x0 = x;
        x_ini = 2*(center*(x'/R)*center/R - x) + x;
        iter_old = 1;
        obj_min = 9999999;
        x_min = [];
        
        while(1)
            x = randn(1,3);
            x = x/norm(x)*R;
            dis = [];
            for iter = 1:M
                dis = [dis norm(XYZ(iter,:) - x)];
            end
            if all(dis<max_dis)&all(dis>min_dis)
                break
            end
        end
        x0 = x;
        beta = zeros(1,M);
        for i =1:1
            for k = 1:M
                beta = solve_eq(F,R,Rb,Rm,Ym,beta,XYZ,x,k);
            end
        end
        beta(beta<0) = 0;
        beta(beta>upper) = upper;
        
        for iter = 1:100000
            [A B C] = ABC(F,R,Rb,Rm,Ym,beta);
            [P, D] = PD(A,B,C,beta,R,Rb);
            obj_old = (G*P'-tau')'*inv_Omega*(G*P'-tau') + delta*ones(M,1)'*P';
            [graP, graD] = graPD(A,B,C,beta,R,Rb);
            %(G*P'-tau')'*inv_Omega*(G*P'-tau')
            S = XYZ;
            dBeta = zeros(M, 2);
            for i = 1:M
                dBeta(i,1) = -(S(i,1) - x(1))/R/sin(D(i)/R)/graD(i);
                dBeta(i,2) = -(S(i,2) - x(2))/R/sin(D(i)/R)/graD(i);
                dBeta(i,3) = -(S(i,3) - x(3))/R/sin(D(i)/R)/graD(i);
            end
%             [hessP, hessD] = hessPD(A,B,C,beta,R,Rb);
%             invH = diag(1./diag(hessP));
%             dP_x = (2*G'*inv_Omega*(G*P'-tau')+ delta*ones(length(P),1))'.*(invH*graP')'*dBeta;
            dP_x = (2*G'*inv_Omega*(G*P'-tau') + delta*ones(length(P),1))'.*graP*dBeta;
            while(1)
                if ss*max(abs(dP_x))>thres
                    ss = 0.5*ss;
                else
                    break
                end
            end
            x = x - ss*dP_x;
            if isnan(x(1))
                break
            end
            ss = ss*alpha;
            for i =1:1
                for k = 1:M
                    beta = solve_eq(F,R,Rb,Rm,Ym,beta,XYZ,x,k);
                end
            end
            if ~all(beta>0)
                if norm(x(1:2))>=R
                    x(1:2) = x(1:2)/norm(x(1:2))*R;
                end
                x(3) = - sqrt(R^2 - (x(1)^2+x(2)^2));
                for i =1:1
                    for k = 1:M
                        beta = solve_eq(F,R,Rb,Rm,Ym,beta,XYZ,x,k);
                    end
                end
            end
            %     beta(beta<0) = 0;
            %     beta(beta>upper) = upper;
            [A B C] = ABC(F,R,Rb,Rm,Ym,beta);
            [P, D] = PD(A,B,C,beta,R,Rb);
            obj = (G*P'-tau')'*inv_Omega*(G*P'-tau') + delta*ones(M,1)'*P';
            if obj<1000
                thres = 1;
                alpha = 1.005;
            end
            dis = norm(x - emitter');
            if plt == 1
                scatter3(x(1),x(2),x(3),5,'m')
            end
            fprintf("obj:%2.2f dis:%2.2f step size:%2.6f\n",obj,dis,ss);
            
            if obj_min>obj
                iter_old = iter;
                obj_min = obj;
                x_min = x;
            end
            if abs(obj)<1e-7||(abs(iter-iter_old)>200)
                break
            end
        end
        x = x_min;
        obj = obj_min;
        dis = norm(x - emitter');
        %Output results
        if plt == 1
            scatter3(emitter(1),emitter(2),emitter(3),50,'filled','r')
            scatter3(x(1),x(2),x(3),50,'k*')
            %text(emitter(1),emitter(2),emitter(3),'e')
        end
        for i = 1:M
            if plt == 1
                scatter3(XYZ(i,1),XYZ(i,2),XYZ(i,3),'filled','c')
                temp = strcat('s ',num2str(i));
                text(XYZ(i,1),XYZ(i,2),XYZ(i,3),temp);
            end
            fprintf("Sensor %d:(%2.2f,%2.2f,%2.2f)\n",i,XYZ(i,1),XYZ(i,2),XYZ(i,3))
        end
        fprintf("Non-Perturbated:r0:%2.2f km,Rm:%2.2f km,Ym:%2.2f km,fc:%2.2f MHz,f:%2.2f MHz\n",R,Rm,Ym,fc,f)
        fprintf("Perturbated: r0:%2.2f km,Rm:%2.2f km,Ym:%2.2f km,fc:%2.2f MHz,f:%2.2f MHz\n",R,P_Rm,P_Ym,P_fc,P_f)
        fprintf("True Location:(%2.2f,%2.2f,%2.2f)\n",emitter(1),emitter(2),emitter(3))
        fprintf("Location:(%2.2f,%2.2f,%2.2f)\n",x(1),x(2),x(3))
        fprintf("True Flying angle:(")
        for i = 1:M-1
            fprintf("%2.2f,",beta0(i))
        end
        fprintf("%2.2f)\n",beta0(M))
        fprintf("Flying angle:(")
        for i = 1:M-1
            fprintf("%2.2f,",beta(i))
        end
        fprintf("%2.2f)\n",beta(M))
        [lg lt] = XYZtoLGLT(x(1),x(2),x(3),R);
        fprintf("The distance to the emitter is: %2.2f km\n",dis)
        fprintf("%s\n",'Localization successful!')
        % [x beta obj] = Alter_Pro_Grad(M,N,P_F,R,P_Rb,P_Rm,P_Ym,G,tau,inv_Omega,delta,upper,max_dis,min_dis,XYZ,plt);
        % diff = norm(emitter'-x,2)
        fprintf(fid,"%2.2f,%d,%2.2f\n",noise_level,index,dis);
    end
end
fprintf(fid,"\t\n");
fclose(fid);