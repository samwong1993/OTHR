%% Projection on a ball
%created by Huang Sen
%Email: huangsen1993@gmail.com
function x = proj_ball(cen,r,x)
    x = r*(x - cen )/norm(x -cen,2)+cen;
end