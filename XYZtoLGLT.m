%% %% (x,y,z) to Latitude and Longitude
%created by Huang Sen
%Email: huangsen1993@gmail.com
function [lon,lat] = XYZtoLGLT(x,y,z,radius)
    lt= asin(z/radius);
    lg = atan(y/x);
    lon = radtodeg(lg);
    lon(lon<0) = 180 + lon;
    lat = radtodeg(lt);
    lat(lat<0) = 180 + lat;
end