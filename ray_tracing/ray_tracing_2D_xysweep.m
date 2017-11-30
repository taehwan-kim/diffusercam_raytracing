header;
addpath ..;
%% fix the random seed
rand('state',2016);
randn('state',2016);

%% parameters:
% size of the diffuser: xrange
% index of the material: np
% number of x/theta bins: numberofxbins
% distance btw diffuser & sensor: z0
% distance btw object & diffuser: z1
% diffuser height std dev.: nstd
% diffuser kernel FWHM: fwhm in the number of pixels

% in Nick's ICCP paper:
% ? = 18?m, with rms surface height of 1.15?m
% pixel pitch: 6.5um
% z0 = 648um
% n=1, np=1.5

lambda = 0.532;
sensorpixel = 6.5;

indexEnv = 1;
indexDiff = 1.5;
z0 = [5000 10000 20000 30000];
% z0 = 5000;
% z1 = (4:0.2:6);
% z1 = fliplr(z1);
% z1 = 1e-5*z1;
% z1 = round(1./z1);
z1 = 20000;
% z1 = [10000 16667 100000];
% z1 = 10000;
% offset = [0 500 1000 1500 2000 2500 3000];
offset = (0:10:150);
% offset = 1000;
nstd = 1.15;
fwhm_raw = 18;

% fwhm = round(fwhm_raw/xpixel);
% sigma = fwhm/2.355;
% sensorrange = sensorpixel * (numberofsensorbins - 1) / 2;
% Fresnel = (sigma*xpixel)^2./(z0*lambda);

%% size of the object in number of voxels
vx = 1;
vy = 1;
vz = 1;

%% the voxel spec in the object space
voxX = 1;
voxY = 1;
voxZ = 100;

%% number of rays from the object
rays = 50000000;
raysPerVoxel = round(rays./(vx*vy*vz));

%% angle spread in degrees
thetaSpread = 20;
phiSpread = 20;

%% number of planes in z-direction for each voxel
k = 1;
dz = voxZ/k;

%% sensor dimensions
npx =1024;
npy = 1024;
x_range = sensorpixel * npx;
y_range = sensorpixel * npy;
center = voxX * vx / 2;
% center = 0;

%% define x & theta grid
% xsense = sensorrange * transpose(linspace(-1, 1, numberofsensorbins));
% xin = xrange * transpose(linspace(-1, 1, numberofxbins));
% thetain = thetarange * transpose(linspace(-1, 1, numberofthetabins));
% 
% sensoredges = cat(1, xsense(1:end) - sensorpixel/2, xsense(end)+sensorpixel/2);

%% get diffuser surface model

% alpha = (fwhm*5-1)/(2*sigma);
% Dheight = nstd*randn(1,length(xin)+1);
% w = gausswin(fwhm*5, alpha);
% y = filter(w,1,Dheight);
% yd = (y(1,2:end)-y(1,1:end-1))/(xpixel);
% ydmat = repmat(yd',[numberofthetabins,1]);
strength = 1;
in = load('./Output/half_deg.mat');
diffuser = in.filtered * strength;
x = in.x; % diffuser coordinate system in physical units (um)
% x = x - min(x);
y = x;
px = mean(diff(x)); % diffuser pixel size in physical units (um)
[Fx, Fy] = gradient(diffuser);

%% define measurement matrix
hMatrix = zeros(npx*npy, vx*vy*vz);
psfStack = zeros(length(offset), length(z1), npx, npy);
psfStack_linear = zeros(length(offset), length(z1), npx*npy);
% xsense = linspace(-1-(sensepixel)/2,1+sensorpixel/2,((length(xin)+1)/subsample));
% xsense = linspace(
% xsense = xrange * transpose(linspace(-1-1/(numberofxbins-1), 1+1/(numberofxbins-1), (0.5*numberofxbins+1)+1));

% caustic = zeros(length(xsense),length(z0));
% 
% xinmat = repmat(xin,[numberofthetabins,1]);
% thetainmat = repelem(thetain, numberofxbins);

%% modeling light field

% plane wave
% inten_unit = zeros(numberofthetabins,1);
% inten_unit((numberofthetabins+1)/2) = 1;
% intensity = repelem(inten_unit,numberofxbins);

% from point source
% intensity = reshape(eye(numberofthetabins),numberofxbins*numberofxbins,1);
% intensity = ones(numberofxbins*numberofxbins,1);

%% refracted by diffuser & propagation to the sensor
% thetaout = (n * thetainmat/np + (n/np-1) * ydmat);
% 
% xout = repmat(xinmat,[1,length(z0)]) + thetaout * z0;
%% accumulate the beam that falls into the sensor pixel

count = 0;
for zzz = 1:length(z0)
    for xxx = 1:length(offset)
        for a = 1:vx
            for b = 1:vy
                for c = 1:vz
                    gatherer = zeros(npy, npx);
                    for j = 1:k
                        % generate random rays from the voxel of interest
                        % everything is in physical units (degrees, ums)
                        xr = (rand(raysPerVoxel,1) + (a-1)) * voxX;
                        yr = (rand(raysPerVoxel,1) + (b-1)) * voxY;
                        th = (rand(raysPerVoxel,1) - 0.5) * thetaSpread;
                        ph = (rand(raysPerVoxel,1) - 0.5) * phiSpread;

                        % distance from z-plane to diffuser
                        z = z1 + voxZ * (c-1) + dz * (j-1);

                        % propagate rays to the diffuser
                        xo = z * tand(th) + xr - center + offset(xxx);
                        yo = z * tand(ph) + yr - center;

                        % get diffuser gradient @ the ray-hitting positions
                        if x(1) < floor(min(xo))-5
                            Fyr = interp2(x(ceil(length(x)/2)+floor(min(xo))-5:ceil(length(x)/2)+ceil(max(xo))+6),y(ceil(length(y)/2)+floor(min(yo))-5:ceil(length(y)/2)+ceil(max(yo))+6),Fy((ceil(length(y)/2)+floor(min(xo))-5:ceil(length(y)/2)+ceil(max(xo))+6),(ceil(length(y)/2)+floor(min(xo))-5:ceil(length(y)/2)+ceil(max(xo))+6)),xo,yo);
                            Fxr = interp2(x(ceil(length(x)/2)+floor(min(xo))-5:ceil(length(x)/2)+ceil(max(xo))+6),y(ceil(length(y)/2)+floor(min(yo))-5:ceil(length(y)/2)+ceil(max(yo))+6),Fx((ceil(length(y)/2)+floor(min(xo))-5:ceil(length(y)/2)+ceil(max(xo))+6),(ceil(length(y)/2)+floor(min(xo))-5:ceil(length(y)/2)+ceil(max(xo))+6)),xo,yo);
                        else
                            Fyr = interp2(x,y,Fy,xo,yo);
                            Fxr = interp2(x,y,Fx,xo,yo);
                        end

                        % throwing out points that did not hit the diffuser and could
                        % not be interpolated
                        good = ~isnan(Fxr) & ~isnan(Fyr);
                        Fxr = Fxr(good);
                        Fyr = Fyr(good);
                        th = th(good);
                        ph = ph(good);
                        xo = xo(good);
                        yo = yo(good);

                        %outputs direction cosines after refraction
                        [uxp, uyp, uzp] = refraction(Fxr, Fyr, th, ph, indexDiff, indexEnv,'angles');

                        %propagate from diffuser to the sensor
                        [yo, xo] = propagation(uyp, uzp, z0(zzz), yo, uxp, xo);

                        %create image at the sensor using 2D histogramming
                        dpx = x_range/npx;
                        dpy = y_range/npy;
                        gatherer = gatherer + hist4(xo,yo,npx,npy,dpx,dpy,x_range/2,y_range/2,px);
                        %imagesc(gatherer);
                        %pause(1/24);
                    end
                    hMatrix(:,(sub2ind([vy,vx,vz],b,a,c))) = gatherer(:);

                    %count to display progress to user
                    count = count + 1;
                    [num2str(count) '/' num2str(vx*vy*vz)];
                end
            end
        end
        psfStack(zzz,xxx,:,:) = reshape(sum(hMatrix,2),npx,npy);
        psfStack_linear(zzz,xxx,:) = reshape(squeeze(psfStack(zzz,xxx,:,:)),1,npx*npy);
    end
end

% anchor = 10;
% autocor = zeros(1,length(z1));
% for i=1:length(z1)
%     temp = corrcoef(psfStack_linear(anchor,:),psfStack_linear(i,:));
%     autocor(i) = temp(1,2);
% end
% plot(autocor);
%figure();
%subplot(3,1,1);
%imagesc(squeeze(psfStack(1,1,:,:)));
%colormap(cm_viridis);
%colorbar;
%subplot(3,1,2);
%imagesc(squeeze(psfStack(1,5,:,:)));
%colormap(cm_viridis);
%colorbar;
%subplot(3,1,3);
%imagesc(squeeze(psfStack(1,10,:,:)));
%colormap(cm_viridis);
%colorbar;
figure(1);
imagesc(squeeze(psfStack(1,1,:,:)),[0,100]);
% colorbar;
colormap(cm_viridis);
axis off;

% anchor = round(length(z1)/2);
anchor = 1;
% anchor = length(z1);
autocor = zeros(length(z0),length(offset));

for j=1:length(z0)
    for i=1:length(offset)
        temp = corrcoef(psfStack_linear(j,anchor,:),psfStack_linear(j,i,:));
        autocor(j,i) = temp(1,2);
    end
end

figure(2);
hold on;
for j=1:length(z0)
    plot(offset,autocor(j,:));
end
% 
% figure();
% imagesc(squeeze(psfStack(4,end,:,:)),[0 10]);
% axis off;
% set(gca,'Position',[0 0 1 1])

% 
% figure(2);
% subplot(3,1,1);
% imagesc(squeeze(psfStack(1,1,:,:)),[0 500]);
% colormap(cm_viridis);
% colorbar;
% subplot(3,1,2);
% imagesc(squeeze(psfStack(1,end/2,:,:)),[0 250]);
% colormap(cm_viridis);
% colorbar;
% subplot(3,1,3);
% imagesc(squeeze(psfStack(1,end,:,:)),[0 20]);
% colormap(cm_viridis);
% colorbar;

% t1 = reshape(squeeze(psfStack(1,:,:)),1,npx*npy);
% t2 = reshape(squeeze(psfStack(2,:,:)),1,npx*npy);
% t3 = reshape(squeeze(psfStack(3,:,:)),1,npx*npy);

% classify = zeros(size(xout,1), size(xout,2));
% 
% for j=1:size(xout,2)
%     classify(:,j) = discretize(xout(:,j),sensoredges);
%     for i=1:size(xout,1)
%         if ~isnan(classify(i,j))
%             caustic(classify(i,j),j) = caustic(classify(i,j),j) + intensity(i);
%         end
%     end
% end
% 
% caustic = caustic + 0.01*randn(size(caustic,1),size(caustic,2));

%% summarize the result
% figure(1)
% subplot(2,1,1);
% hold on;
% plot(Dheight);
% plot(y);
% legend('dh','y');
% subplot(2,1,2);
% plot(yd);
% legend boxoff;
% figure(2)
% for i=1:length(z0)
%     subplot(length(z0),1,i);
%     plot(xsense,caustic(:,i))
%     xlim([-1000 1000]);
%     ylim([-1 10]);
%     title(strcat('F=',string(Fresnel(i))));
% end

% xsense = zeros(1,length(xin));

