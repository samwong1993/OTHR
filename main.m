clear
R = 6371.2;
%plt = 1 plot the earth,emitter,sensors and the generated sequence
plt = 1;
if plt == 1
    figure('color','k')
    ha=axes('units','normalized','position',[0 0 1 1]);
    uistack(ha,'down')
    img=imread('background.jpg');
    image(img)
    colormap gray
    set(ha,'handlevisibility','off','visible','off');
    rotate3d on
    hold on
    %Add legend
    %emitter
    point1 = scatter3(0,0,0,50,'filled','r');
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
% beta0 = [0.349025176895742,0.279243732945529,0.174632834143026,0.183395025597251,0.285385790728298];  
% XYZ = zeros(M,3);  
% [x0 y0 z0] = LGLTtoXYZ(116.24,39.55,R);  
% emitter = [x0 y0 z0]';  
% [x0 y0 z0] = LGLTtoXYZ(128.72,40.55,R);  
% XYZ(1,:) = [x0 y0 z0];  
% [x0 y0 z0] = LGLTtoXYZ(130.42,38.68,R);  
% XYZ(2,:) = [x0 y0 z0];  
% [x0 y0 z0] = LGLTtoXYZ(132.94,33.82,R);  
% XYZ(3,:) = [x0 y0 z0];  
% [x0 y0 z0] = LGLTtoXYZ(130.90,31.84,R);  
% XYZ(4,:) = [x0 y0 z0];  
% [x0 y0 z0] = LGLTtoXYZ(129.06,35.63,R);  
% XYZ(5,:) = [x0 y0 z0];  


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
tau = generate_tau(M,F,R,Rb,Rm,Ym,emitter,XYZ);
%Step 3
[x beta obj] = GPGD(M,N,P_F,R,P_Rb,P_Rm,P_Ym,G(1:M-1,:),tau(1:M-1),inv_Omega,upper,max_dis,min_dis,XYZ,plt);
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
%Lengend
if plt == 1
	h = legend([point1,point2,point3,point4,point5],'Emitter', 'Sensors', 'Initial Point','Estimated Location','Generated Points','AutoUpdate','off');
    %set(h,'box','off')
end
%plot ray path
if plt == 1
    for i = 1:M
        raypath(emitter,XYZ(i,:),Rm,Rb,Ym,R);
    end
end