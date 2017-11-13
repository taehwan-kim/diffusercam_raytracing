header;

% parameters:
% size of the diffuser: xrange
% index of the material: np
% number of x/theta bins: numberofxbins
% distance btw diffuser & sensor: z0
% distance btw object & diffuser: z1
% diffuser height std dev.: nstd
% diffuser kernel FWHM: fwhm pixels

n = 1;
xrange = 1;
np = 1.5;
numberofxbins = 31;
z0 = 0.5;
z1 = 8;
nstd = 0.01;
fwhm = 10;

thetarange = xrange / z1;
numberofthetabins = numberofxbins;

xin = xrange * transpose(linspace(-1, 1, numberofxbins));
yin = xrange * transpose(linspace(-1, 1, numberofxbins));

Dheight = nstd*randn(length(xin)+1,length(xin)+1);
Dfilter = imgaussfilt(Dheight,fwhm);
% w = gausswin(fwhm);
% y = filter(w,1,Dheight);
Gy = (Dfilter(:,2:end)-Dfilter(:,1:end-1))/(xin(2)-xin(1));
Gx = (Dfilter(2:end,:)-Dfilter(1:end-1,:))/(xin(2)-xin(1));
% yd = (y(1,2:end)-y(1,1:end-1))/(xin(2)-xin(1));
% ydmat = repmat(yd',[numberofthetabins,1]);

xpixel = xin(2)-xin(1);
subsample = 10;
sensepixel = xpixel * subsample;
xsense = linspace(-1-(sensepixel)/2,1+sensepixel/2,((length(xin)+1)/subsample));
ysense = linspace(-1-(sensepixel)/2,1+sensepixel/2,((length(xin)+1)/subsample));
% [XS,YS] = meshgrid(xsense,ysense);
% xsense = xrange * transpose(linspace(-1-1/(numberofxbins-1), 1+1/(numberofxbins-1), (0.5*numberofxbins+1)+1));

% caustic = zeros(1,length(xsense)-1);
thetain = thetarange * transpose(linspace(-1, 1, numberofthetabins));
phiin = thetarange * transpose(linspace(-1, 1, numberofthetabins));
caustic = zeros(length(xsense)-1,length(xsense)-1);

% hist(thetain*z1,xin);

xinmat = repmat(xin,[numberofthetabins,1]);
thetainmat = repelem(thetain, numberofxbins);
yinmat = repmat(yin,[numberofthetabins,1]);
phiinmat = repelem(phiin, numberofxbins);

intensityx_vec = reshape(eye(numberofthetabins),numberofxbins*numberofxbins,1);
intensityy_vec = reshape(eye(numberofthetabins),numberofxbins*numberofxbins,1);
intensityx = repmat(intensityx_vec,[1,size(intensityy_vec)]);
intensityy = repmat(intensityy_vec',[size(intensityx_vec),1]);
% intensity = ones(numberofxbins*numberofxbins,1);


yout(1,:) = repmat(xinmat,[1,numberofxbins]) + z0 * (n * thetainmat/np + (n/np-1) * Dx(Gx(),xinmat));
plot(xout);
classify = discretize(xout,xsense);

for i=1:length(xout)
    if ~isnan(classify(i))
        caustic(classify(i)) = caustic(classify(i)) + intensity(i);
    end
end

figure(2)
plot(xsense(1:end-1)+sensepixel/2,caustic)
% xsense = zeros(1,length(xin));


function y = Dx(ydmat,xin)
    y = ydmat.*xin;
end