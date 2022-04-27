% script to extract grounding line depths from BedMachinev4
% for Slater 2022 GRL
clear; close all;

% bedmachine dataset is v4 available from https://nsidc.org/data/IDBMG4

% load bedmachine dataset
bmfile = '~/Documents/BedMachineGreenland-2021-04-20.nc';
x = single(ncread(bmfile,'x'));
y = single(ncread(bmfile,'y'));
t = single(ncread(bmfile,'mask'))';
b = single(ncread(bmfile,'bed'))';

% subsample if wanting to speed up code (ss=1 means no subsampling)
ss = 1;
x = x(1:ss:end);
y = y(1:ss:end);
t = t(1:ss:end,1:ss:end);
b = b(1:ss:end,1:ss:end);
io = 0*b;
t2 = t;
[X,Y] = meshgrid(x,y);

% extract ice-ocean boundary pixels
for i=2:size(t,1)-1,
    for j=2:size(t,2)-1,
        if t(i,j)==2,
            kern = t(i-1:i+1,j-1:j+1);
            if sum(kern(:)==0 | kern(:)==3)>0,
                io(i,j)=1;
                t2(i,j)=5;
            end
        end
    end
end

% in a few cases the bedmachine mask merges two calving fronts
% that have to be manually separated
io(find(y==-909125),find(x==-2225)) = 0;
io(find(y==-1587725),find(x==-361925)) = 0;
io(find(y==-1339025),find(x==545575)) = 0;
io(find(y==-2295575),find(x==491125)) = 0;
io(find(y==-2297075),find(x==493975)) = 0;
io(find(y==-2684975),find(x==228925)) = 0;

% group pixels into calving fronts (i.e. connected pixels)
calvingfronts = bwconncomp(io);

% place results into convenient structure
for i=1:length(calvingfronts.PixelIdxList),
    cfs(i).inds = calvingfronts.PixelIdxList{i};
    cfs(i).bed = b(cfs(i).inds);
    [cfs(i).gldepth,glind] = min(cfs(i).bed);
    cfs(i).x = X(cfs(i).inds(glind));
    cfs(i).y = Y(cfs(i).inds(glind));
end

% save
save calvingfronts.mat cfs

% optional sense check plot
figure(); hold on;
imagesc(x,y,t2); axis xy equal;
plot([cfs.x],[cfs.y],'ko');
for i=1:length([cfs.x]),
    text(cfs(i).x,cfs(i).y,num2str(i),'fontsize',16);
end



