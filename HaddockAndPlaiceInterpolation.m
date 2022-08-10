function HaddockAndPlaiceInterpolation
% loading landmark data and declaring landmark variables
img1 = imread("fish01edit.jpg");
img2 = imread("fish02edit.jpg");

load("fin1data.mat");
load("fin2data.mat");
load("fin3data.mat");
load("fin4data.mat");
load("fin5data.mat");
load("fishmouthdata.mat")
load("fisheyedata.mat")
xi = [fin1x;fin2x;fin3x;fin4x;fin5x];%;mouthx;eyex];
yi = [fin1y;fin2y;fin3y;fin4y;fin5y];%;mouthy;eyey];

Xi = [fin1X;fin2X;fin3X;fin4X;fin5X];%;mouthX;eyeX];
Yi = [fin1Y;fin2Y;fin3Y;fin4Y;fin5Y];%;mouthY;eyeY];

fineness = 5; % grid spacing

% ----GENERATING MATRIX EQAUTIONS----
m = length(xi); % number of points mapped
A = zeros(m);

for i = 1:m 
    for j = 1:m
        A(i,j) = gk(xi(i,1),yi(i,1),xi(j,1),yi(j,1));
    end
end

F1 = linsolve(A,Xi(:,1));
H1 = linsolve(A,Yi(:,1));

[coord1,coord2] = meshgrid(1:fineness:352,1:fineness:274);

% ----CALCULATING SURFACES----
u_surf = 0;
for k = 1:m
    u_surf = u_surf + F1(k)*gk(coord1,coord2,xi(k,1),yi(k,1));
end

v_surf = 0;
for k = 1:m
    v_surf = v_surf + H1(k)*gk(coord1,coord2,xi(k,1),yi(k,1));
end

% plotting mappings as surfaces in their own right
surf(coord1,coord2,u_surf)
hold on
for i = 1:m
    plot3(xi(i,1),yi(i,1),Xi(i,1),"o","MarkerFaceColor","k",'MarkerEdgeColor','k',"MarkerSize",5)
end
hold off
title('x coords')

surf(coord1,coord2,v_surf)
hold on
for i = 1:m
    plot3(xi(i,1),yi(i,1),Yi(i,1),"o","MarkerFaceColor","k",'MarkerEdgeColor','k',"MarkerSize",5)
end
hold off
title('y coords')

% ----BENDING ENERGY----
bending_energy_u = transpose(F1(1:m))*A*F1(1:m)
bending_energy_v = transpose(H1(1:m))*A*H1(1:m)

% ----CALCULATING MAPPINGS----
B = ones(m,3);

for i = 1:m
    B(i,2) = xi(i,1);
    B(i,3) = yi(i,1);
end

C = [A,B;transpose(B),zeros(3)];

F = linsolve(C,[Xi(1:m,1);0;0;0]); 
H = linsolve(C,[Yi(1:m,1);0;0;0]);

u = 0;
for k = 1:m
    u = u + F(k)*gk(coord1,coord2,xi(k),yi(k));
end


% adding affine part
u = u + F(m+1) + F(m+2)*coord1 + F(m+3)*coord2;

v = 0;
for k = 1:m
    v = v + H(k)*gk(coord1,coord2,xi(k),yi(k));
end

v = v + H(m+1) + H(m+2)*coord1 + H(m+3)*coord2;

% ----CALCULATING LINE ELEMENTS----
uprimex = 0;
for k = 1:m
    uprimex = uprimex + F(k)*gkprimex(coord1,coord2,xi(k),yi(k));
end
uprimex = uprimex + F(m+2);

uprimey = 0;
for k = 1:m
    uprimey = uprimey + F(k)*gkprimey(coord1,coord2,xi(k),yi(k));
end
uprimey = uprimey + F(m+3);

vprimex = 0;
for k = 1:m
    vprimex = vprimex + H(k)*gkprimex(coord1,coord2,xi(k),yi(k));
end
vprimex = vprimex + H(m+2);

vprimey = 0;
for k = 1:m
    vprimey = vprimey + H(k)*gkprimey(coord1,coord2,xi(k),yi(k));
end
vprimey = vprimey + H(m+3);

%JM = [uprimex,uprimey;vprimex,vprimey]; % jacobian matrix

J = uprimex.*vprimey - uprimey.*vprimex; % the jacobian determinant

dX = ones(size(coord1)) * fineness; % defining how long dX and dY are
dY = ones(size(coord2)) * fineness;

dx = (uprimex .* dX) + (uprimey .* dY); % applying jacobian component wise to all dX and dY
dy = (vprimex .* dX) + (vprimey .* dY);

% ----UNTRANSFORMED PLOTS----
% Plotting untransformed fish and grid
warp(coord1,coord2,zeros(size(u)),img1) 
hold on
surf(coord1,coord2,zeros(size(u)),'FaceAlpha',0)
hold off

% plotting landmarks
hold on
for k = 1:m
    plot3(xi(k,1),yi(k,1),0,"x","MarkerSize",15,"MarkerEdgeColor",'r',"LineWidth",2);
end
hold off

%visuals
view(0,90)
set(gca,'FontSize',20)
pbaspect([1,1,1])
set(gca,'YDir','normal');
set(get(gca,'ylabel'),'rotation',0)
axis equal tight
xlabel('x')
ylabel('y')

imshow(img2)
% plotting landmarks
hold on
for k = 1:m
    plot3(Xi(k,1),Yi(k,1),0,"x","MarkerSize",15,"MarkerEdgeColor",'b',"LineWidth",2);
end
hold off

%visuals
view(0,90)
set(gca,'FontSize',20)
pbaspect([1,1,1])
set(gca,'YDir','normal');
set(get(gca,'ylabel'),'rotation',0)
axis equal tight

% ----TRANSFORMED PLOTS----
% Plotting transformed fish and grid
warp(u,v,zeros(size(u)),img1) 
hold on
surf(u,v,zeros(size(u)),'FaceAlpha',0)
hold off

% plotting landmarks
hold on
for k = 1:m
    plot3(Xi(k,1),Yi(k,1),0,"x","MarkerSize",15,"MarkerEdgeColor",'b',"LineWidth",2);
end
hold off

% plotting line elements
hold on
for k = 1:m
    %quiver(u,v,dX,dY,0.5,"Color",'b')
end
hold off

%visuals
view(0,90)
set(gca,'FontSize',20)
pbaspect([1,1,1])
set(gca,'YDir','normal');
set(get(gca,'ylabel'),'rotation',0)
axis equal tight
xlabel('x')
ylabel('y')

%----D'ARCY THOMPSON STYLE SUPERIMPOSITIONS----
% desired fish with image grid
imshow(img2)
hold on
surf(u,v,zeros(size(u)),'FaceAlpha',0)
hold off

%visuals
view(0,90)
set(gca,'FontSize',20)
pbaspect([1,1,1])
set(gca,'YDir','normal');
set(get(gca,'ylabel'),'rotation',0)
axis equal tight

%----DEFORMATION ANALYSIS----
pcolor([1:fineness:352],[1:fineness:274],J)
shading interp;
colorbar;
hold on
for k = 1:m
    plot3(xi(k,1),yi(k,1),0,"x","MarkerSize",15,"MarkerEdgeColor",'r',"LineWidth",2);
end
hold off

%----SIDE BY SIDE----
warp(u,v,zeros(size(u)),img1) 
%visuals
view(0,90)
set(gca,'FontSize',20)
pbaspect([1,1,1])
set(gca,'YDir','normal');
set(get(gca,'ylabel'),'rotation',0)
axis equal tight
set(gca,'visible','off')

imshow(img2)
%visuals
view(0,90)
set(gca,'FontSize',20)
pbaspect([1,1,1])
set(gca,'YDir','normal');
set(get(gca,'ylabel'),'rotation',0)
axis equal tight
end

