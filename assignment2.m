clear
clc

nx = 20;
ny = 30;
Vo = 2;
delx = 1e-3;
dely = 1e-3;

G = sparse(nx*ny,nx*ny);
Vv = zeros(nx*ny,1);
V = zeros(nx,ny);
B = zeros((nx*ny),1);


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
%             ryp = (cMap(i.j) + cMap(i,j+1))/2;            
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
%             ryp = (cMap(i.j) + cMap(i,j+1))/2;            
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
%             ryp = (cMap(i.j) + cMap(i,j+1))/2;            
%             rxm = (cMap(i,j) + cMap(i-1,j))/2;        
%             rxp = (cMap(i,j) + cMap(i+1,j))/2;       
      
            G(n,n) = -4;%-(rxm+rxp+rym+ryp);
            G(n,nym) = 1;%rym;
            G(n,nyp) = 1;%ryp;
            G(n,nxm) = 1;%rxm;
            G(n,nxp) = 1;%rxp; 
            B(n,1) = 0;
        end
   
    end
end

spy(G)
Vv = G\B;

for j = 1:ny
    for i = 1:nx
        n = i+(j-1)*nx;
        V(i,j) = Vv(n,1);
    end
end
        

surf(V)


























