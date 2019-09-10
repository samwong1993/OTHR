%% Generate G
%created by Huang Sen
%Email: huangsen1993@gmail.com
function [G] = generate_G(N,M)
%Generate G
G = zeros(N,M);
m = 1;
k = 2;
for i = 1:N
    G(i,m) = 1;
    G(i,k) = -1;
    k = k + 1;
    if k > M
        m = m + 1;
        k = m + 1;
    end
end
end