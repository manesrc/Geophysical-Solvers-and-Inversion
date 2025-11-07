function [X_units,T] = CreaMalla_rectangulo(nen, L_ref, InfoMesh)
% [X,T] = CreateMesh(RefElement,dom,nx,ny)
% Creates the topology of an structured and uniform mesh over a rectangular
% domain
%
% Input:
%   RefElement: reference element
%   nen:  number of element nodes
%   dom = [x1,x2,y1,y2]:  vertices' coordinates
%   nx,ny: number of elements in each direction
% Output:
%   X:  nodal coordinates
%   T:  connectivities

elem = InfoMesh.elemType; % (only option implemented)
x1 = InfoMesh.ini_x;
x2 = InfoMesh.fin_x;
y1 = InfoMesh.ini_y;
y2 = InfoMesh.fin_y;
nx = InfoMesh.nel_x; 
ny = InfoMesh.nel_y;

if elem == 1 % quads
    if nen == 9
        degree = 2;
    else 
        degree =1;
    end
else
    sprintf('element not implemented')
end

npx = nx*degree+1;
npy = ny*degree+1;
npt = npx * npy;
x = linspace(x1,x2,npx); 
y = linspace(y1,y2,npy); 
[x,y] = meshgrid(x,y); 
X = [reshape(x',npt,1), reshape(y',npt,1)]; 
X_units = [X(:,1) X(:,2)] * L_ref;

% Connectivities
if elem == 1            % Quadrilaterals
    if nen == 4         % Q1
        T = zeros(nx*ny,4);
        for i=1:ny
            for j=1:nx
                ielem = (i-1)*nx+j;
                inode = (i-1)*(npx)+j;
                T(ielem,:) = [inode   inode+1   inode+npx+1   inode+npx];
            end
        end
    elseif nen == 9     % Q2
        T = zeros(nx*ny,9);
        for i=1:ny
            for j=1:nx
                ielem = (i-1)*nx + j;
                inode = (i-1)*2*npx + 2*(j-1) + 1;
                T(ielem,:)=[inode   inode+2   inode+2*(npx)+2   inode+2*(npx) ...
                    inode+1   inode+2+(npx)   inode+2*(npx)+1   inode+(npx)   inode+(npx)+1 ];
            end
        end
    end
elseif elem == 2        % Triangles
    if nen == 3             % P1
        nx_2 = nx/2; ny_2 = ny/2; 
        
        T = zeros(2*nx*ny,3);
        for i=1:ny
            for j=1:nx
                ielem = 2*((i-1)*nx+j)-1;
                inode = (i-1)*(npx)+j;
                nodes = [inode   inode+1   inode+npx+1    inode+npx]; 
                if (i<=ny_2 && j<=nx_2) || (i>ny_2 && j>nx_2)
                    T(ielem,:) = nodes([1,2,3]); 
                    T(ielem+1,:) = nodes([1,3,4]); 
                else
                    T(ielem,:) = nodes([1,2,4]); 
                    T(ielem+1,:) = nodes([2,3,4]); 
                end
            end   
        end
    elseif nen == 6         % P2
        nx_2 = round(nx/2); ny_2 = round(ny/2);         
        T = zeros(2*nx*ny,6); 
        for i=1:ny
            for j=1:nx
                ielem=2*((i-1)*nx+j)-1;
                inode=(i-1)*2*(npx)+2*(j-1)+1;
                nodes = [inode+(0:2)  inode+npx+(0:2)  inode+2*npx+(0:2)];
                if (i<=ny_2 && j<=nx_2) || (i>ny_2 && j>nx_2)
                    T(ielem,:)   = nodes([1  9  7  5  8  4]); 
                    T(ielem+1,:) = nodes([1  3  9  2  6  5]); 
                else
                    T(ielem,:)   = nodes([1  3  7  2  5  4]);
                    T(ielem+1,:) = nodes([3  9  7  6  8  5]);
                end
            end    
        end
    end

end   


end 