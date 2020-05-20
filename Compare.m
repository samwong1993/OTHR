fid = 'randM5WTtau0_0001.txt';
[noise_level,index,diff]=textread(fid,'%f%d%f','delimiter',',');
for i = 1:length(index)
        error(index(i),round(noise_level(i))) =  diff(i);
end
error(error<20)=0;
error(error>20)=1;
excel = sum(error);
zz = 100 -excel;
y=[zz;excel]';
b=bar(y);
grid on;
ch = get(b,'children');
set(gca,'XTickLabel',{'0','100','200','300','400','500','600','700','800','900','1000'})
legend('Localization Success','Localization Failed');
legend boxoff;
% xlabel('Standard Deviation of TDOA Measurement Noise (ns)')
% ylabel('RMSE of the Geolocation(m)')
xlabel('Standard Deviation of TDOA Measurement Noise (ns)')
ylabel('Number');