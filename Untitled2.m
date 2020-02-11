clear
clc

W = 20;
L = 30;
Vo = 5;
del = 1e-3;

G = sparse(W,L);
Vcalc = zeros((W*L),1);
B = zeros((W*L),1);


%Initialize Boundary Conditions
for i = 1:W
    for j = 1
        B(i,1) = Vo;
    end
end
        
 
for j = 1:L
    for i = 1:W
        G(i,j) = n;
        G(i-1,j) = nxm;
        G(i+1,j) = nxp;
    end
end

 V = B\G




























