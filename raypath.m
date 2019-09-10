%% Plot ray path
%created by Huang Sen
%Email: huangsen1993@gmail.com
function raypath(emitter,S,Rm,Rb,Ym,R)
    t = 0:0.1:1;
    dis = norm(emitter'-S)/2;
    b = dis;
    a=b/dis^2;
    R_max =  Rm*Rb/(Rb-Ym);
    x0 = (t-0.5)*dis*2;
    y0 = -a*x0.^2+b;
    for i =1:length(t)
        xyz(i,:) = emitter' + t(i)*(S-emitter');
        xyz(i,:) = xyz(i,:)/norm(xyz(i,:))*(R+y0(i));
    end
    comet3(xyz(:,1),xyz(:,2),xyz(:,3),0.1);
end

