clear
clc

nx = 20;
ny = 30;
Vo = 1;
delx = 1e-3;
dely = 1e-3;

G = sparse(nx*ny,nx*ny);
Vv = zeros(nx*ny,1);
V = zeros(nx,ny);
B = zeros((nx*ny),1);
        
%Part 1 a)        
        
%Initialize Left Boundary Conditions
for i = 1:nx*ny
    B(i,1) = Vo;
    B(i,ny) = 0;
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
            %Part a) B(n,1) = 0
            B(n,1) = 1;
            
        elseif i == 1
            %Mapping
            nym = (i)+(j-2)*nx;
            nyp = (i)+(j)*nx;
            nxp = (i+1)+(j-1)*nx;
            G(n,n) = 1;
            B(n,1) = 0;

            
        elseif i == nx
            %Mapping
            nym = (i)+(j-2)*nx;
            nyp = (i)+(j)*nx;
            nxm = (i-1)+(j-1)*nx;
            G(n,n) = 1;
            B(n,1) = 0;
                      

        else
            %Mapping
            nym = (i)+(j-2)*nx;
            nyp = (i)+(j)*nx;
            nxm = (i-1)+(j-1)*nx;
            nxp = (i+1)+(j-1)*nx;



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
xlabel('ny Value')
ylabel('nx Value')




%Part 1 b)
% %Analytical Soln
nx = 20;
ny = 30;
W = 3;
L = 2;
Vanal = zeros(ny,nx);
y = linspace(0,L,ny);
x = linspace(-W/2,W/2, nx); 
a = L;
b = W/2;
[X,Y] = meshgrid(x,y);

for n = 1:2:100
    Vanal = Vanal +((4*Vo)/pi)* ((1/n)*(cosh((n*pi*X)./a)./(cosh((n*pi*b)./a))).*sin(((n*pi*Y)./a)));
end


figure(3)
surf(Vanal)
xlabel('ny Value')
ylabel('nx Value')




















%Part 2
nx = 20;
ny = 30;

cMap = zeros(nx, ny);

for k = 1:10
Acond = 1;
Bcond = linspace(10E-10,1,10);
Bcond = Bcond(k);

for j = 1:ny
    for i = 1:nx
        cMap(i,j) = Acond;
        if ((j<((2/3)*ny))&&((j>((1/3)*ny))&&(i<((1/3)*nx))))|| ((j<((2/3)*ny))&&((j>((1/3)*ny))&&(i>((2/3)*nx))))
            cMap(i,j) = Bcond;
        end
    end
end
figure(4)
surf(cMap)    
xlabel('ny Value')
ylabel('nx Value')

%Initialize Left Boundary Conditions
for i = 1:nx*ny
    B(i,1) = Vo;
    B(i,ny) = 0;
end
        
%Set diagonal
for j = 1:ny
    for i = 1:nx
        n = i+(j-1)*nx;
        %Set the Boundary Nodes
        if j == 1
            G(n,:) = 0;
            G(n,n) = 1; 
%             B(n,1) = 0;
            
        elseif j == ny
            G(n,:) = 0;
            G(n,n) = 1;
            B(n,1) = 0;
            
        elseif i == 1
            %Mapping
            nxm = (i)+(j-2)*nx;
            nxp = (i)+(j)*nx;
            nyp = (i+1)+(j-1)*nx;
            G(n,:) = 0;
            G(n,n) = 1;
            B(n,1) = 0;
            
            rxm = ((cMap(i,j) + cMap(i,j-1))/2);
            rxp = ((cMap(i,j) + cMap(i,j+1))/2);            
            ryp = ((cMap(i,j) + cMap(i+1,j))/2);

            G(n,n) = -(rxp+rxm+ryp);
            G(n,nxm) = rxm;
            G(n,nyp) = ryp;
            G(n,nxp) = rxp;
            
        elseif i == nx
            %Mapping
            nym = (i-1)+(j-1)*nx;
            nxp = (i)+(j)*nx;
            nxm = (i)+(j-2)*nx;
            G(n,:) = 0;
            G(n,n) = 1;
            B(n,1) = 0;
            
            rxm = ((cMap(i,j) + cMap(i,j-1))/2);
            rxp = ((cMap(i,j) + cMap(i,j+1))/2);            
            rym = ((cMap(i,j) + cMap(i-1,j))/2);

            G(n,n) = -(rxm+rym+rxp);
            G(n,nym) = rym;
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;            

        else
            %Mapping
            nym = (i-1)+(j-1)*nx;
            nyp = (i+1)+(j-1)*nx;
            nxm = (i)+(j-2)*nx;
            nxp = (i)+(j)*nx;

            rym = ((cMap(i,j) + cMap(i-1,j))/2);
            ryp = ((cMap(i,j) + cMap(i+1,j))/2);            
            rxm = ((cMap(i,j) + cMap(i,j-1))/2);        
            rxp = ((cMap(i,j) + cMap(i,j+1))/2);       
      
   
            G(n,n) = -(rxm+rxp+rym+ryp);
            G(n,nym) = rym;
            G(n,nyp) = ryp;
            G(n,nxm) = rxm;
            G(n,nxp) = rxp; 
            B(n,1) = 0;
            


        end
   
    end
end

Vv = G\B;

for j = 1:ny
    for i = 1:nx
        n = i+(j-1)*nx;
        V(i,j) = Vv(n,1);
    end
end

figure(5)
surf(V)
xlabel('ny Value')
ylabel('nx Value')



%Gradient of V
for j = 1:ny
    for i = 1:nx
        if j == 1
            Ex(i,j) = (V(i,j+1)-V(i,j));
        elseif j == ny
            Ex(i,j) = (V(i,j)-V(i,j-1));
        else
            Ex(i,j) = (V(i,j+1)-V(i,j-1))*0.5;
        end
        if i == 1
            Ey(i,j) = (V(i+1,j)-V(i,j));
        elseif i == nx
            Ey(i,j) = (V(i,j)-V(i-1,j));
        else
            Ey(i,j) = (V(i+1,j)-V(i-1,j))*0.5;
        end
    end
end

Ex = -Ex;
Ey = -Ey;
% figure(6)
% surf(Ex)
% xlabel('ny Value')
% ylabel('nx Value')
% 
% figure(7)
% surf(Ey)
% xlabel('ny Value')
% ylabel('nx Value')
% 
% figure(8)
% surf(Ex,Ey)
% xlabel('ny Value')
% ylabel('nx Value')

Jx = cMap.*Ex;
Jy = cMap.*Ey;

% figure(9)
% surf(Jx)
% xlabel('ny Value')
% ylabel('nx Value')
% 
% figure(10)
% surf(Jy)
% xlabel('ny Value')
% ylabel('nx Value')
% 
% figure(11)
% surf(Jx,Jy)
% xlabel('ny Value')
% ylabel('nx Value')

C0 = sum(Jx(1,:));
Cnx = sum(Jx(nx,:));
Curr = (C0+Cnx)*0.5;

figure(15)
plot(Bcond, Curr,'r.')
xlabel('Bcond (sigma)')
ylabel('Current (A)')
hold on
end



















