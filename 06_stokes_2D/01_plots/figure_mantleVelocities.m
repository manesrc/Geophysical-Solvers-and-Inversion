function figure_mantleVelocities(X,velocities,Interphase,fig_num)
% plot the velocity results results 
% INPUTS: 
% X: Mesh
% velocities: mantle flow results
% Interphase: definition of the interphase
% fig_num: number of the figure

n_nodes = size(X,1);
velocities = reshape(velocities,[n_nodes 2]);

V_x = velocities(:,1);
V_y = velocities(:,2);

%% figure
figure(fig_num);
clf
hold on
Xmin = min(X(:,1)); Xmax = max(X(:,1));
[cx,cy] = mySlice(Xmin,Xmax,grad,n1,nz);
conepl1 = coneplot(X(:,1),X(:,2),V_x,V_y,cx,cy);





function [cx,cy,cz] = mySlice(Xmin,Xmax,grad,n1,nz)
%
% solo funciona para gradientes con 2 componentes negativos... falta
% agregar los otros cuadrantes... seguramente se puede hacer menos ifs con
% un poco de gracia...
%
if nargin < 5
   nz = n1;
end
xmin = Xmin(1);
xmax = Xmax(1);
xmid = (xmax-xmin)/2;

ymin = Xmin(2);
ymax = Xmax(2);
ymid = (ymax-ymin)/2;

zmin = Xmin(3);
zmax = Xmax(3);


gradx = grad(1);
grady = grad(2);



if grad(1)<=0 && grad(2)<=0
   if grad(1) > grad(2)
      tanAlpha = gradx / grady;
      
      x1 = xmid + tanAlpha*ymid;
      y1 = ymax;
      x2 = xmid - tanAlpha*ymid;
      y2 = ymin;
   else
      tanAlpha = grady / gradx;
      x1 = xmax;
      y1 = ymid + tanAlpha*xmid;
      x2 = xmin;
      y2 = ymid - tanAlpha*xmid;
   end
elseif grad(1)>=0 && grad(2)<=0
   if grad(1) > grad(2)
      tanAlpha = gradx / grady;
      
      x1 = xmid + tanAlpha*ymid;
      y1 = ymax;
      x2 = xmid - tanAlpha*ymid;
      y2 = ymin;
   else
      tanAlpha = grady / gradx;
      x1 = xmax;
      y1 = ymid - tanAlpha*xmid;
      x2 = xmin;
      y2 = ymid + tanAlpha*xmid;
   end
elseif grad(1)<=0 && grad(2)>=0
   % faltan los otros cuadrantes
else  % grad(1)>=0 && grad(2)>=0
   % faltan los otros cuadrantes
end

xlin = linspace(x1,x2,n1);
ylin = linspace(y1,y2,n1);
zlin = linspace(zmin,zmax,nz);

cx = zeros(n1,nz);
cy = cx;
cz = cx;
for I = 1:n1
   for J = 1:nz
      cx(I,J) = xlin(I); 
      cy(I,J) = ylin(I); 
      cz(I,J) = zlin(J); 
   end
end

end

end
