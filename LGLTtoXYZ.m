%% Latitude and Longitude to (x,y,z)
%created by Huang Sen
%Email: huangsen1993@gmail.com
function [x y z] = LGLTtoXYZ(lon,lat,radius)
    lg = degtorad(lon);
    lt = degtorad(lat);
    temp = radius*cos(lt);
    x = temp*cos(lg);
    y = temp*sin(lg);
    z = radius*sin(lt);
end