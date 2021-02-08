clc
clear
close all

for i = 1:10
    for j = 1:10
        z(i,j) = i^2+j^2;
    end
end
for i =1:10
    for j = 1:10
        X1(10*i+j,1) = i;
        X1(10*i+j,2) = j;
    end
end
X = X1(11:110,:);
X = (-1).*X.^2+10.*X;
Yt = reshape(z,100,1);
Wout = pinv(X'*X+(1e-8).*eye(size(X,2)))*(X'*Yt);

for i = 1:10
    for j = 1:10
        T(i,j) = ((-1).*[i,j].^2+10.*[i,j]) * Wout;
    end
end

figure
mesh(T)