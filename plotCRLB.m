%% Plot RMSE 
%created by Huang Sen
%Email: huangsen1993@gmail.com
clf
clear
%%plot monte carlo results for CRLB
figure(1)
pltMCCRLB('GDM3_100.txt','*k-')
hold on
pltMCCRLB('GDM5_100.txt','^r-')
pltMCCRLB('ICASSPM5_0_1.txt','og-')
%%plot CRLB
R = 6371.2;
Rm = 6650;
Ym = 100;
Rb = Rm - Ym;
fc = 10;
f = 15;
F = f/fc;
%Input location of emitter and sensors
beta0 = [0.114957231412252,0.449398124172348,0.277420425918117,0.0168095219080640,0.103488345084960];
%Hong Kong
[x0 y0 z0] = LGLTtoXYZ(114.16,22.28,R);
emitter = [x0 y0  z0]';
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

%add noise(ns)
sigma = [0:100:1000];
c = 3*10^5;
sigma_axis = sigma;
%sigma = sigma_axis;
sigma = 10^(-9)*sigma*c;
crlb = zeros(length(sigma),1);
M = 3;
index = [3 4 5] ;
for noise_level = 2:length(sigma)
    crlb(noise_level) = CRLB_tdoaOTHR(F, Rb, Ym, Rm, R, beta0, XYZ, emitter, sigma(noise_level), M,index);
end
plot(sigma_axis,1000*crlb,'-', 'linewidth', 1.1, 'markerfacecolor', [29, 191, 151]/255);

M = 5;
index = [1 2 3 4 5] ;
for noise_level = 2:length(sigma)
    crlb(noise_level) = CRLB_tdoaOTHR(F, Rb, Ym, Rm, R, beta0, XYZ, emitter, sigma(noise_level), M,index);
end
plot(sigma_axis,1000*crlb,'r-', 'linewidth', 1.1, 'markerfacecolor', [36, 169, 225]/255);


grid on
legend('GPGD Monte Carlo simulation with 3 sensors','GPGD Monte Carlo simulation with 5 sensors','Quasi-Newton Monte Carlo simulation with 5 sensors','CRLB with 3 sensors','CRLB with 5 sensors')