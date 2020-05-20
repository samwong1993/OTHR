function x_min = Initialization(XYZ,max_dis,min_dis,M,inv_Omega,G,tau,F,R,Rb,Rm,Ym,upper,x_input)  
while(1)
    x = randn(1,3);
    x = x/norm(x)*R;
    dis = [];
    for iter = 1:M
        dis = [dis norm(XYZ(iter,:) - x)];
    end
    if all(dis<max_dis)&all(dis>min_dis)
        break
    end
end
x0 = x;
beta = zeros(1,M);
for i =1:1
    for k = 1:M
        beta = solve_eq(F,R,Rb,Rm,Ym,beta,XYZ,x,k);
    end
end
beta(beta<0) = 0;
beta(beta>upper) = upper;
[A B C] = ABC(F,R,Rb,Rm,Ym,beta);
[P, D] = PD(A,B,C,beta,R,Rb);
obj_x0 = (G*P'-tau')'*inv_Omega*(G*P'-tau');
center = R*sum(XYZ)/M/(norm(sum(XYZ)/M));
x_ini = 2*(center*(x'/R)*center/R - x) + x;
for i =1:1
    for k = 1:M
        beta = solve_eq(F,R,Rb,Rm,Ym,beta,XYZ,x_ini,k);
    end
end
beta(beta<0) = 0;
beta(beta>upper) = upper;
[A B C] = ABC(F,R,Rb,Rm,Ym,beta);
[P, D] = PD(A,B,C,beta,R,Rb);
obj_ini = (G*P'-tau')'*inv_Omega*(G*P'-tau');
if obj_ini>obj_x0
    obj_min = obj_x0;
    x_min = x0;
else
    obj_min = obj_ini;
    x_min = x_ini;
end

if ~isempty(x_input)
    beta = zeros(1,M);
    for i =1:1
        for k = 1:M
            beta = solve_eq(F,R,Rb,Rm,Ym,beta,XYZ,x_input,k);
        end
    end
    beta(beta<0) = 0;
    beta(beta>upper) = upper;
    [A B C] = ABC(F,R,Rb,Rm,Ym,beta);
    [P, D] = PD(A,B,C,beta,R,Rb);
    obj_input = (G*P'-tau')'*inv_Omega*(G*P'-tau');
    if obj_min>obj_input
        obj_min = obj_input;
        x_min = x_input;
    end
end

end