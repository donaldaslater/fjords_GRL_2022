% script to define boundary between fjords and shelf for Slater 2022 GRL
clear; close all;

% bedmachine dataset is v4 available from https://nsidc.org/data/IDBMG4

%% load bedmachine dataset

bmfile = '~/Documents/BedMachineGreenland-2021-04-20.nc';
x = single(ncread(bmfile,'x'));
y = single(ncread(bmfile,'y'));
t = single(ncread(bmfile,'mask'))';
b = single(ncread(bmfile,'bed'))';
% make t into a mask in which 0 is water and 1 is greenland land/ice
t(t>=3) = 0;
t(t>1) = 1;
% subsample dataset to make code faster (ss=1 means no subsampling)
ss = 5;
x = x(1:ss:end);
y = y(1:ss:end);
t = t(1:ss:end,1:ss:end);
b = b(1:ss:end,1:ss:end);

% due to subsample, can get bits of water that are unconnected to the ocean
% we want to set these to be land/ice
not_t = 1-t;
a = bwconncomp(not_t,8);
for i=1:length(a.PixelIdxList),
    l(i) = length(a.PixelIdxList{i});
end
[~,id] = max(l);
inds = [];
for i=1:length(a.PixelIdxList),
    if l(i)~=l(id),
        inds = [inds;a.PixelIdxList{i}];
    end
end
t(inds) = 1;

%% define fjord/shelf boundary and masks

% remove islands
% calculate connected land area (thus isolating islands)
a = bwconncomp(t,8);
for i=1:length(a.PixelIdxList),
    l(i) = length(a.PixelIdxList{i});
end
% get mainland
[~,id] = max(l);
% get (x,y) coords of mainland pixels
inds = a.PixelIdxList{id};
[X,Y] = meshgrid(x,y);
xb = X(inds); yb = Y(inds);

% do shrinkwrap to separate fjord and shelf
% s=0.8 seems good for ss=10, s=0.5 for ss=5
s = 0.5;
k = boundary(double(xb),double(yb),s);
xb = xb(k); yb = yb(k);

% % optional plotting to check sensitivity to tightness parameter
% s = [0.25,0.5,0.75]; % tightness parameter
% figure(); hold on;
% imagesc(x,y,t); axis xy  equal;
% for i=1:length(s),
%     k = boundary(double(xb),double(yb),s(i));
%     xplot = xb(k); yplot = yb(k);
%     p(i) = plot(xplot,yplot,'linewidth',1);
% end
% legend(p,{'1','2','3'});

% make fjords mask f
% mask is 0 outside of perimeter (i.e. shelf/ocean)
% mask is 0.5 inside of perimeter if water (i.e. fjords)
% mask is 1 inside of perimeter if not water (i.e. land/ice)
f = t;

% those that are water
inds1 = find(t==0);

% of which are coastal (use manual polygon to speed up code)
xc = 1e5*[0.5784,1.5970,1.9070,2.7485,3.1914,4.7415,5.6715,6.4244,...
          7.8416,8.5502,8.9931,8.0188,8.3731,6.8673,5.9372,5.9815,...
          3.8114,1.8185,-0.3516,-2.4332,-4.1604,-4.9576,-6.2862,...
          -6.9062,-6.6848,-6.0648,-5.0019,-4.4704,-4.3375,-4.2490,...
          -4.1604,-4.4704,-4.4704,-4.5147,-4.2047,-3.8504,-3.2746,...
          -2.6546,-1.3703];
yc = 1e6*[-3.3704,-3.3305,-3.1135,-2.9452,-2.7637,-2.6618,-2.4669,...
          -2.3695,-2.2632,-2.1126,-1.9975,-1.7185,-1.4439,-1.1914,...
          -0.9567,-0.7397,-0.6556,-0.6423,-0.6998,-0.7884,-0.8991,...
          -1.0054,-1.0984,-1.2047,-1.3819,-1.4882,-1.4837,-1.5546,...
          -1.7318,-1.9266,-2.0595,-2.2233,-2.3784,-2.5644,-2.7504,...
          -2.9187,-3.0781,-3.2464,-3.3305];
mi = inpolygon(X(inds1),Y(inds1),xc,yc);
inds1(mi==0) = [];

% of which lie inside boundary
mi = inpolygon(X(inds1),Y(inds1),xb,yb);
inds1(mi==0) = [];
f(inds1) = 0.5;

% binary fjords mask fb (fjords=1 and everything else=0)
fb = f;
fb(find(f==1))=0;
fb(find(f==0.5))=1;

% ice mask m
m = single(ncread(bmfile,'mask'))';
m = m(1:ss:end,1:ss:end);
m(find(m~=2)) = NaN;

% % optional plot
% figure(); hold on;
% imagesc(x,y,f); axis xy equal;
% plot(xb,yb,'k');

%% get fjord effective depths

shelfinds = find(f==0); % linear indices of shelf points
ed = NaN(size(b)); % initialise effective depth

% cycle through z
z = [0,-25,-50:-50:-2000];
for j=1:length(z),
    
    % everything shallower than z(j) is assigned 0
    % everything that is deeper than or equal to z(j) is assigned 1
    b0 = 0*b;
    inds = find(b<=z(j));
    b0(inds) = 1;

    % get connectivity map of points assigned with 1
    a = bwconncomp(b0);

    % for each connected group, check if is connected to shelf
    % if is connected, set effective depth to z(j) because there is a
    % continuous path at this depth from the edge to these points
    % by cycling through z from the surface downwards we end up with the
    % deepest connected depth
    for i=1:a.NumObjects,
        if ~isempty(find(ismember(a.PixelIdxList{i},shelfinds))),
            ed(a.PixelIdxList{i}) = z(j);
        end
    end
    
end

ed = single(ed);

%% save

save fjords.mat x y f fb t b xb yb m ed




