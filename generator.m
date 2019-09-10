%% Generate emitter and sensors
%created by Huang Sen
%Email: huangsen1993@gmail.com
function [emitter,XYZ,beta0] = generator(M,F,R,Rb,Rm,Ym,max_dis,min_dis)
emitter = randn(3,1);
emitter = R*emitter / norm(emitter);
for i = 1:M
while(1)
    XYZ(i,:) = randn(1,3);
    XYZ(i,:) = R*XYZ(i,:) / norm(XYZ(i,:));
    if norm(emitter' - XYZ(i,:))<max_dis&norm(emitter' - XYZ(i,:))>min_dis
        break
    end
end
end
%Calculate corresponing beta of x
while(1)
beta = zeros(1,M);
x = emitter';
for i =1:10
    for k = 1:M
        beta = solve_eq(F,R,Rb,Rm,Ym,beta,XYZ,x,k);
    end
end
[A B C] = ABC(F,R,Rb,Rm,Ym,beta);
[P D] = PD(A,B,C,beta,R,Rb);
for k = 1:M
    penalty(k) = norm(XYZ(k,:)-x,2)^2 - 4*R^2*(1 - cos(D(k)/R))/2;
end
beta0 = beta;
if sum(penalty) < 1e-5
    break;
end
end
end