function KiteMappings
m = 4; % number of points mapped
A = zeros(m);
fineness = 0.2;
% starting coords
X1 = [0,1]; % top 
X2 = [1,0]; % right
X3 = [-1,0]; % left
X4 = [0,-1]; % bottom

% mapped to coords
x1 = [0,1]; % top
x2 = [1,0]; % right
x3 = [-1,0]; % left
x4 = [0,-2]; % bottom

Xi = [X1;X2;X3;X4];
xi = [x1;x2;x3;x4];

% generating matrix of eqns
for i = 1:m 
    for j = 1:m
        A(i,j) = gk(Xi(i,1),Xi(i,2),Xi(j,1),Xi(j,2));
    end
end

[coord1,coord2] = meshgrid(-3:fineness:3,-3:fineness:3);
[coord1surf,coord2surf] = meshgrid(-5:fineness:5,-5:fineness:5);

F1 = linsolve(A,xi(:,1));
H1 = linsolve(A,xi(:,2));

% ----CALCULATING SURFACES----
u_surf = 0;
for k = 1:m
    u_surf = u_surf + F1(k)*gk(coord1surf,coord2surf,Xi(k,1),Xi(k,2));
end

v_surf = 0;
for k = 1:m
    v_surf = v_surf + H1(k)*gk(coord1surf,coord2surf,Xi(k,1),Xi(k,2));
end

% plotting surfaces
figure
subplot(1,2,1)
surf(coord1surf,coord2surf,u_surf)
hold on
for i = 1:m
    plot3(Xi(i,1),Xi(i,2),xi(i,1),"o","MarkerFaceColor","k",'MarkerEdgeColor','k',"MarkerSize",10)
end
hold off
title('x coordinate surface')
xlabel('x')
ylabel('y')
set(gca,'fontsize',15)


subplot(1,2,2)
surf(coord1surf,coord2surf,v_surf)
hold on
for i = 1:m
    plot3(Xi(i,1),Xi(i,2),xi(i,2),"o","MarkerFaceColor","k",'MarkerEdgeColor','k',"MarkerSize",10)
end
hold off
title('y coordinate surface')
xlabel('x')
ylabel('y')
set(gca,'fontsize',15)
% ----CALCULATING MAPPINGS----
B = ones(m,3);

for i = 1:m
    B(i,2) = Xi(i,1);
    B(i,3) = Xi(i,2);
end

C = [A,B;transpose(B),zeros(3)];

F = linsolve(C,[xi(1:m,1);0;0;0]); % xi(1:m) x coords
H = linsolve(C,[xi(1:m,2);0;0;0]); % xi(1:m) y coords

u = 0;
for k = 1:m
    u = u + F(k)*gk(coord1,coord2,Xi(k,1),Xi(k,2));
end

% adding affine part
u = u + F(m+1) + F(m+2)*coord1 + F(m+3)*coord2;

v = 0;
for k = 1:m
    v = v + H(k)*gk(coord1,coord2,Xi(k,1),Xi(k,2));
end

v = v + H(m+1) + H(m+2)*coord1 + H(m+3)*coord2;

% ----BENDING ENERGY----
bending_energy_u = transpose(F(1:m))*A*F(1:m)
bending_energy_v = transpose(H(1:m))*A*H(1:m)


% ----CALCULATING LINE ELEMENTS----
uprimex = 0;
for k = 1:m
    uprimex = uprimex + F(k)*gkprimex(coord1,coord2,Xi(k,1),Xi(k,2));
end
uprimex = uprimex + F(m+2);

uprimey = 0;
for k = 1:m
    uprimey = uprimey + F(k)*gkprimey(coord1,coord2,Xi(k,1),Xi(k,2));
end
uprimey = uprimey + F(m+3);

vprimex = 0;
for k = 1:m
    vprimex = vprimex + H(k)*gkprimex(coord1,coord2,Xi(k,1),Xi(k,2));
end
vprimex = vprimex + H(m+2);

vprimey = 0;
for k = 1:m
    vprimey = vprimey + H(k)*gkprimey(coord1,coord2,Xi(k,1),Xi(k,2));
end
vprimey = vprimey + H(m+3);

% defining how long dX and dY are
dX = ones(size(coord1)) * fineness;
dY = zeros(size(coord2)) * fineness;

% line element transformations
dx = (uprimex .* dX) + (uprimey .* dY);
dy = (vprimex .* dX) + (vprimey .* dY);

% the jacobian determinant
J = uprimex.*vprimey - uprimey.*vprimex; 

dx_stretch = zeros(size(dX));
dy_stretch = zeros(size(dY));
dx_rotate = zeros(size(dX));
dy_rotate = zeros(size(dY));
x_principal_stretch = zeros(size(dX));
y_principal_stretch = zeros(size(dY));
x1_principal_stretch_vec = zeros(size(dX));
y1_principal_stretch_vec = zeros(size(dY));
x2_principal_stretch_vec = zeros(size(dX));
y2_principal_stretch_vec = zeros(size(dY));
volume_change = zeros(size(dX));
% Jacobian matrix at a certain coordinate (i.e. certain area element)
for i = 1:size(dX)
    for j = 1:size(dY)
        JM = [uprimex(i,j),uprimey(i,j);vprimex(i,j),vprimey(i,j)]; % Jacobian matrix
        c = transpose(JM)*JM; % cauchy-green deformation gradient tensor
        [eigvecs,eigvals] = eig(c);
        x_principal_stretch(i,j) = sqrt(eigvals(1,1));
        x1_principal_stretch_vec(i,j) = eigvecs(1,1);
        y1_principal_stretch_vec(i,j) = eigvecs(1,2);
        x2_principal_stretch_vec(i,j) = eigvecs(2,1);
        y1_principal_stretch_vec(i,j) = eigvecs(2,2);
        y_principal_stretch(i,j) = sqrt(eigvals(2,2));
        volume_change(i,j) = x_principal_stretch(i,j)*y_principal_stretch(i,j); % volume change using princp stretches
        U = sqrtm(c); % stretch matrix
        R = U\JM; % rotation matrix
        dx_stretch(i,j) = U(1,1)*dX(i,j) + U(1,2)*dY(i,j); % stretching part of line elements
        dy_stretch(i,j) = U(2,1)*dX(i,j) + U(2,2)*dY(i,j);
        dx_rotate(i,j) = R(1,1)*dX(i,j) + R(1,2)*dY(i,j); % rotation part of line elements
        dy_rotate(i,j) = R(2,1)*dX(i,j) + R(2,2)*dY(i,j);
    end
end


% ----UNTRANSFORMED PLOTS----
% plotting kite
figure
surf(coord1,coord2,zeros(size(u)),'FaceAlpha',0.01)
grid off
hold on
for k = 1:m
    plot3(Xi(k,1),Xi(k,2),0,"o","MarkerFaceColor","k",'MarkerEdgeColor','k',"MarkerSize",5);
end
% plotting line elements
%quiver(coord1,coord2,dX,dY,1,"Color",'b')
hold off

% visuals and axis
view(0,90)
set(gca,'FontSize',20)
pbaspect([1,1,1])
xlim([-3,3])
ylim([-3,3])


% ----TRANSFORMED PLOTS----
% plotting transformed kite
figure
surf(u,v,zeros(size(u)),'FaceAlpha',0.01)
hold on
for k = 1:m
    plot3(xi(k,1),xi(k,2),0,"o","MarkerFaceColor",'k','MarkerEdgeColor',"k","MarkerSize",5);
end
% plotting line element changes
%quiver(u,v,dX,dY,'Color','k') % undeformed line element
%quiver(u,v,dx,dy,"Color",'r') % plotting overall deformed line element change
quiver(u,v,dx_rotate,dy_rotate,"Color",'b') % plotting deformed line element rotations
%quiver(u,v,dx_stretch,dy_stretch,'Color','g') % plotting deformed line element stretches

% plotting principal stretches
%quiver(u,v,x1_principal_stretch_vec,y1_principal_stretch_vec,"Color",'b')
%quiver(u,v,x2_principal_stretch_vec,y2_principal_stretch_vec,"Color",'r')

hold off

% visuals and axis
view(0,90)
set(gca,'FontSize',20)
pbaspect([1,1,1])
xlabel('x')
ylabel('y')
set(get(gca,'ylabel'),'rotation',0)
grid off
xlim([-3,3])
ylim([-3,3])

% plotting volume change via jacobian
figure
pcolor(-3:fineness:3,-3:fineness:3,J);
hold on
%quiver(u,v,dx,dy,1,"Color",'r')
for k = 1:m
    plot3(Xi(k,1),Xi(k,2),0,"o","MarkerFaceColor","k",'MarkerEdgeColor','k',"MarkerSize",5);
end
hold off
shading interp;
colorbar
grid off
pbaspect([1,1,1])
xlabel('x')
ylabel('y')
set(get(gca,'ylabel'),'rotation',0)
set(gca,'FontSize',20)

% plotting volume change via principal stretches (indentical to the jacobian of course)
figure
pcolor(-3:fineness:3,-3:fineness:3,volume_change);
hold on
for k = 1:m
    plot3(Xi(k,1),Xi(k,2),0,"o","MarkerFaceColor","k",'MarkerEdgeColor','k',"MarkerSize",5);
end
hold off
shading interp;
colorbar
grid off
pbaspect([1,1,1])
xlabel('x')
ylabel('y')
set(get(gca,'ylabel'),'rotation',0)
set(gca,'FontSize',20)

end

