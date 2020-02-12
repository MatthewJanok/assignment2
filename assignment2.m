clear
clc

nx = 20;
ny = 30;
Vo = 2;
delx = 1e-3;
dely = 1e-3;
Acond = 1;
Bcond = 10E-2;

G = sparse(nx*ny,nx*ny);
Vv = zeros(nx*ny,1);
V = zeros(nx,ny);
B = zeros((nx*ny),1);

cMap = zeros(nx*ny, nx*ny);

for j = 1:ny
    for i = 1:nx
        cMap(i,j) = Acond;
        if ((ny<15)&&(ny>5)&&(nx>20))|| ((ny<15)&&(ny>5)&&(nx<10))
            cMap(i,j) = Bcond;
        end
    end
end
        
        
        
        
%Initialize Left Boundary Conditions
for i = 1:nx
    B(i,1) = Vo;

end
        
%Set diagonal
for j = 1:ny
    for i = 1:nx
        n = i+(j-1)*nx;
        %Set the Boundary Nodes
        if j == 1
            G(n,:) = 0;
            G(n,n) = 1; 
            B(n,1) = 1;
            
        elseif j == ny
            G(n,:) = 0;
            G(n,n) = 1;
            B(n,1) = 1;
            
        elseif i == 1
            %Mapping
            nym = (i)+(j-2)*nx;
            nyp = (i)+(j)*nx;
            nxp = (i+1)+(j-1)*nx;
            G(n,n) = 1;
            B(n,1) = 0;
%             rym = (cMap(i,j) + cMap(i,j-1))/2;
%             ryp = (cMap(i,j) + cMap(i,j+1))/2;            
%             rxp = (cMap(i,j) + cMap(i+1,j))/2;
% 
%             G(n,n) = -(rxp+rym+ryp);
%             G(n,nym) = rym;
%             G(n,nyp) = ryp;
%             G(n,nxp) = rxp;
            
        elseif i == nx
            %Mapping
            nym = (i)+(j-2)*nx;
            nyp = (i)+(j)*nx;
            nxm = (i-1)+(j-1)*nx;
            G(n,n) = 1;
            B(n,1) = 0;
%             rym = (cMap(i,j) + cMap(i,j-1))/2;
%             ryp = (cMap(i,j) + cMap(i,j+1))/2;            
%             rxm = (cMap(i,j) + cMap(i-1,j))/2;
% 
%             G(n,n) = -(rxm+rym+ryp);
%             G(n,nym) = rym;
%             G(n,nyp) = ryp;
%             G(n,nxp) = rxp;            

        else
            %Mapping
            nym = (i)+(j-2)*nx;
            nyp = (i)+(j)*nx;
            nxm = (i-1)+(j-1)*nx;
            nxp = (i+1)+(j-1)*nx;

%             rym = (cMap(i,j) + cMap(i,j-1))/2;
%             ryp = (cMap(i,j) + cMap(i,j+1))/2;            
%             rxm = (cMap(i,j) + cMap(i-1,j))/2;        
%             rxp = (cMap(i,j) + cMap(i+1,j))/2;       
%       
%             
%             G(n,n) = -(rxm+rxp+rym+ryp);
%             G(n,nym) = rym;
%             G(n,nyp) = ryp;
%             G(n,nxm) = rxm;
%             G(n,nxp) = rxp; 
%             B(n,1) = 0;

            G(n,n) = -4;
            G(n,nym) = 1;
            G(n,nyp) = 1;
            G(n,nxm) = 1;
            G(n,nxp) = 1;
            B(n,1) = 0;
        end
   
    end
end

figure(1)
spy(G)
Vv = G\B;

for j = 1:ny
    for i = 1:nx
        n = i+(j-1)*nx;
        V(i,j) = Vv(n,1);
    end
end

figure(2)
surf(V)


%Analytical Soln
Vanal = zeros(nx,ny);
a = linspace(0,ny);
b = linspace(-nx/2,nx/2); 

for n = 1:2:100
    Vanal = Vanal +((4*Vo)/pi)* sum((1/n)*(cosh((n*pi*nx)./a)./(cosh((n*pi*b)./a)).*sin((n*pi*ny./a))));
end

figure(3)
surf(Vanal)


% for j = 1:ny
%     for i = 1:nx
%         if j == 1
%             Ex(i,j) = (V(i,j+1)-V(i,j));
%         elseif j == ny
%             Ex(i,j) = (V(i,j)-V(i,j-1));
%         else
%             Ex(i,j) = (V(i,j+1)-V(i,j-1))*0.5;
%         end
%         if i == 1
%             Ey(i,j) = (V(i+1,j)-V(i,j));
%         elseif i == nx
%             Ey(i,j) = (V(i,j)-V(i-1,j));
%         else
%             Ey(i,j) = (V(i+1,j)-V(i-1,j))*0.5;
%         end
%     end
% end

% Ex = -Ex;
% Ey = -Ex;
% 
% eflowx = cMap.*Ex;
% eflowy = cMap *Ey;
% 
% C0 = sum(eflowx(1,:));
% Cnx = sum(eFLow(nx,:));
% Curr = (C0+Cnx)*0.5;





















