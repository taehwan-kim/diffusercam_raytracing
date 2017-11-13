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
sensorpixel = 1;

indexEnv = 1;
indexDiff = 1.5;
z0 = 1000;
z1 = 50000;
nstd = 1.15;
fwhm_raw = 18;

% fwhm = round(fwhm_raw/xpixel);
% sigma = fwhm/2.355;
% sensorrange = sensorpixel * (numberofsensorbins - 1) / 2;
% Fresnel = (sigma*xpixel)^2./(z0*lambda);

%% size of the object in number of voxels
vx = 10;
vy = 10;
vz = 1;

%% the voxel spec in the object space
voxX = 10;
voxY = 10;
voxZ = 100;

%% number of rays from the object
rays = 50000000;
raysPerVoxel = round(rays./(vx*vy*vz));

%% angle spread in degrees
thetaSpread = 1;
phiSpread = 1;

%% number of planes in z-direction for each voxel
k = 1;
dz = voxZ/k;

%% sensor dimensions
npx = 128;
npy = 128;
x_range = sensorpixel * npx;
y_range = sensorpixel * npy;
center = sensorpixel * npx / 2 - voxX * vx / 2;

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
in = load('../../Output/half_deg.mat');
diffuser = in.filtered * strength;
x = in.x; % diffuser coordinate system in physical units (um)
x = x - min(x);
y = x;
px = mean(diff(x)); % diffuser pixel size in physical units (um)
[Fx, Fy] = gradient(diffuser);

%% define measurement matrix
hMatrix = zeros(npx*npy, vx*vy*vz);
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
                xo = z * tand(th) + xr + center;
                yo = z * tand(ph) + yr + center;
                
                % get diffuser gradient @ the ray-hitting positions
                Fyr = interp2(x(1:ceil(max(xo))+5),y(1:ceil(max(yo))+5),Fy((1:ceil(max(yo))+5),(1:ceil(max(xo))+5)),xo,yo);
                Fxr = interp2(x(1:ceil(max(xo))+5),y(1:ceil(max(yo))+5),Fx((1:ceil(max(yo))+5),(1:ceil(max(xo))+5)),xo,yo);
                
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
                [yo, xo] = propagation(uyp, uzp, z0, yo, uxp, xo);
                
                %create image at the sensor using 2D histogramming
                dpx = x_range/npx;
                dpy = y_range/npy;
                gatherer = gatherer + hist4(xo,yo,npx,npy,dpx,dpy,0,0,px);
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

imagesc(reshape(sum(hMatrix,2),npx,npy));
colormap(cm_viridis);
colorbar;
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

