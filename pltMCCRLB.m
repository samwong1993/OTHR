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
% max_err = [];
% min_err = [];
% for i = 1:n
%     max_err(i) = max(error(:,i));
%     min_err(i) = min(error(:,i));
% end
% z = 10^3*[min_err' max_err'];
% h = area(sigma,z);
% h(1).FaceColor = [1 1 1];
% h(2).FaceColor = [216 233 241]/255;
% set(h,'FaceAlpha',0.5)
% set(h,'EdgeAlpha',0)
% hold on
% pointer = plot(sigma,10^3*err_mean,'ok-','linewidth',1.1,'markerfacecolor',[36, 169, 225]/255);
% legend(pointer,'Mean Value');
% for i = 1:m
%     plot(sigma,10^3*error(i,:),'b--o') 
% end
plot(sigma,10^3*err_mean,fig,'linewidth',1.1);
% xlabel('Standard Deviation of TDOA Measurement Noise (ns)')
% ylabel('RMSE of the Geolocation(m)')
end