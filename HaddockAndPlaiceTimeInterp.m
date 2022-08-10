function  HaddockAndPlaiceTimeInterp
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
xi = [fin1x;fin2x;fin3x;fin4x;fin5x;mouthx;eyex];
yi = [fin1y;fin2y;fin3y;fin4y;fin5y;mouthy;eyey];

oldXi = [fin1X;fin2X;fin3X;fin4X;fin5X;mouthX;eyeX];
oldYi = [fin1Y;fin2Y;fin3Y;fin4Y;fin5Y;mouthY;eyeY];

interplength = 6; % how many sets of intermediate landmarks
interpy = zeros(length(xi),interplength);
interpx = zeros(length(xi),interplength);

%interpolate x coords of landmarks, along each row is a set of x
%coordinates, the columns going down give the next set of x coordinates
%and so on. Here y coordinates are functions of x coordinates
for i = 1:length(xi)
    interpy(i,:) = interp1([xi(i),oldXi(i)],[yi(i),oldYi(i)],linspace(xi(i),oldXi(i),interplength));
end

%interpolate y coords of landmarks in the same way as the x's. We imagine
%the y coordinates as functions of x coordinates.
for i = 1:length(xi)
    interpx(i,:) = interp1([yi(i),oldYi(i)],[xi(i),oldXi(i)],linspace(yi(i),oldYi(i),interplength));
end

fineness = 10;
figure
for i = 1:interplength
    Xi = interpx(:,i);
    Yi = interpy(:,i);
    [u,v] = TPS_NO_LINE_ELEMENTS(xi,yi,Xi,Yi,fineness,img1);
    subplot(2,3,i)
    surf(u,v,zeros(size(u)),'FaceAlpha',0)
    %visuals
    view(0,90)
    set(gca,'FontSize',20)
    pbaspect([1,1,1])
    set(gca,'YDir','normal');
    set(get(gca,'ylabel'),'rotation',0)
    axis equal tight
    set(gca,'visible','off')
end
figure
for i = 1:interplength
    Xi = interpx(:,i);
    Yi = interpy(:,i);
    [u,v] = TPS_NO_LINE_ELEMENTS(xi,yi,Xi,Yi,fineness,img1);
    subplot(2,3,i)
    warp(u,v,zeros(size(u)),img1)
    %visuals
    view(0,90)
    set(gca,'FontSize',20)
    pbaspect([1,1,1])
    set(gca,'YDir','normal');
    set(get(gca,'ylabel'),'rotation',0)
    axis equal tight
    set(gca,'visible','off')
end
    
end

