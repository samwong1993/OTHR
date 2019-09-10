function pltMCCRLB(fid,fig)
[noise_level,index,diff]=textread(fid,'%f%d%f','delimiter',',');
for i = 1:length(index)
        error(index(i),round(noise_level(i))) =  diff(i);
end
sigma = [0:100:1000];
[m n] = size(error);
for i = 1 : n
    err_mean(i) = sqrt(sum(error(:,i).^2)/m);
end
err_mean = mean(error);
plot(sigma,10^3*err_mean,fig, 'linewidth', 1.1)
% for i = 1:m
%     plot(sigma,10^3*error(i,:),'b--o') 
% end
xlabel('Standard Deviation of TDOA Measurement Noise (ns)')
ylabel('RMSE of the Geolocation(m)')
end