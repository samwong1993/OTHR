%% Generalized Projected Gradient Descent Algorithm
%created by Huang Sen
%Email: huangsen1993@gmail.com
%mo is the momentum parameter and distance
function [x beta obj] = GPGD(M,N,F,R,Rb,Rm,Ym,G,tau,inv_Omega,upper,max_dis,min_dis,XYZ,plt)
        thres = 10;
        mo = 1.05;
        ss = 1e-4;
        x = [];
        x = Initialization(XYZ,max_dis,min_dis,M,inv_Omega,G,tau,F,R,Rb,Rm,Ym,upper,x);
        if plt == 1
            scatter3(x(1),x(2),x(3),50,'filled','b')
        end
        beta = zeros(1,M);
        for i =1:20
            for k = 1:M
                beta = solve_eq(F,R,Rb,Rm,Ym,beta,XYZ,x,k);
            end
        end
        iter_old = 1;
        obj_min = 9999999;
        x_min = [];
        for  iter = 1:1000000
            x_old = x;
            if mod(iter,20)==0&&obj>10000
                x = Initialization(XYZ,max_dis,min_dis,M,inv_Omega,G,tau,F,R,Rb,Rm,Ym,upper,x);  
            end
            [A B C] = ABC(F,R,Rb,Rm,Ym,beta);
            [P, D] = PD(A,B,C,beta,R,Rb);
            obj_old = (G*P'-tau')'*inv_Omega*(G*P'-tau');
            [graP, graD] = graPD(A,B,C,beta,R,Rb);
            S = XYZ;
            dBeta = zeros(M, 2);
            for i = 1:M
                dBeta(i,1) = -(S(i,1) - S(i,3)*x(1)/x(3))/R/sin(D(i)/R)/graD(i);
                dBeta(i,2) = -(S(i,2) - S(i,3)*x(2)/x(3))/R/sin(D(i)/R)/graD(i);
            end
            %Newton method
%             [hessP, hessD] = hessPD(A,B,C,beta,R,Rb);
%             invH = diag(1./diag(hessP));
%             dP_x = (2*G'*inv_Omega*(G*P'-tau'))'.*(invH*graP')'*dBeta;
            %Gradient descent
            dP_x = (2*G'*inv_Omega*(G*P'-tau'))'.*graP*dBeta;
            dP_x(3) = 0;
            while(1)
                if ss*max(abs(dP_x))>thres
                    ss = 0.5*ss;
                else
                    break
                end
            end
            x = x - ss*dP_x;
            if isnan(x(1))
                break
            end
            ss = ss*mo;
            if norm(x(1:2))>=R
                x(x>R) = R;
                x(x<-R) = -R;
                x(1:2) = x(1:2)/norm(x(1:2))*R*0.99;
            end
            x(3) = sqrt(R^2 - (x(1)^2+x(2)^2));
            for i =1:20
                for k = 1:M
                    beta = solve_eq(F,R,Rb,Rm,Ym,beta,XYZ,x,k);
                end
            end
            if ~all(beta>0)
                if norm(x(1:2))>=R
                    x(1:2) = x(1:2)/norm(x(1:2))*R;
                end
                x(3) = - sqrt(R^2 - (x(1)^2+x(2)^2));
                for i =1:20
                    for k = 1:M
                        beta = solve_eq(F,R,Rb,Rm,Ym,beta,XYZ,x,k);
                    end
                end
            end
            beta(beta<0) = 0;
            beta(beta>upper) = upper;
            [A B C] = ABC(F,R,Rb,Rm,Ym,beta);
            [P, D] = PD(A,B,C,beta,R,Rb);
            obj = (G*P'-tau')'*inv_Omega*(G*P'-tau');
            if obj<20
                thres = 1;
                mo = 1.005;
            end
            if plt == 1
                scatter3(x(1),x(2),x(3),5,'m')
            end
            fprintf("obj:%2.2f step size:%2.6f\n",obj,ss);
            
            if obj_min>obj
                iter_old = iter;
                obj_min = obj;
                x_min = x;
            end
            if abs(obj)<1e-7||(abs(iter-iter_old)>200&obj<10000)
                break
            end
        end
        x = x_min;
        obj = obj_min;
        fprintf("obj:%2.2f step size:%2.6f\n",obj,ss);
end