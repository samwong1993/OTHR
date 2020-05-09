%% Hildreth's Algorithm
%created by Huang Sen
%Email: huangsen1993@gmail.com
%Projection onto polyhedron
function x = hildreth(XYZ,max_dis,min_dis,M,R,x_0)  
b = zeros(M,1);
for i = 1:M
    A(2*i-1,:) = XYZ(i,:);
    b(2*i-1) = R^2 - min_dis^2/2;
    A(2*i,:) = - XYZ(i,:);
    b(2*i) = max_dis^2/2 - R^2;
end
x = x_0;
z = zeros(2*M,1);
k = 0;
index = [1:2*M];
x_old = randn(3,1);
while(1)
    k = k + 1;
    index_k = mod(k,2*M)+1;
    c = min(z(index_k),0.1*(b(index_k)-A(index_k,:)*x)/norm(A(index_k,:))^2);
    x = x + c*A(index_k,:)';
    z(index_k==i) = z(index_k==i) - c;
    if mod(k,2*M+1)==0
       if x_old == x
           break
       else
           x_old = x;
       end
    end
end
% cvx_begin
%     variable x(3)
%     minimize( norm(x-x_0) )
%     subject to
%         R^2 - max_dis/2<= XYZ(1,:)*x<= R^2 - min_dis/2
%         R^2 - max_dis/2<= XYZ(2,:)*x<= R^2 - min_dis/2
%         R^2 - max_dis/2<= XYZ(3,:)*x<= R^2 - min_dis/2
% cvx_end
end
