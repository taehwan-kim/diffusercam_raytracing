header;

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
n = 1;
np = 1.5;
sensorpixel = 6.5;
xpixel = sensorpixel/2;
numberofsensorbins = 1025;
numberofxbins = 513;
numberofthetabins = 5;
z0 = [10 50 100 200 400 800];
z1 = 10;
nstd = 1.15;
fwhm_raw = 18;

fwhm = round(fwhm_raw/xpixel);
sigma = fwhm/2.355;
sensorrange = sensorpixel * (numberofsensorbins - 1) / 2;
xrange = xpixel * (numberofxbins - 1) / 2;
thetarange = xrange / z1;
Fresnel = (sigma*xpixel)^2./(z0*lambda);

%% define x & theta grid
xsense = sensorrange * transpose(linspace(-1, 1, numberofsensorbins));
xin = xrange * transpose(linspace(-1, 1, numberofxbins));
thetain = thetarange * transpose(linspace(-1, 1, numberofthetabins));

sensoredges = cat(1, xsense(1:end) - sensorpixel/2, xsense(end)+sensorpixel/2);

%% get diffuser surface model

alpha = (fwhm*5-1)/(2*sigma);
Dheight = nstd*randn(1,length(xin)+1);
w = gausswin(fwhm*5, alpha);
y = filter(w,1,Dheight);
yd = (y(1,2:end)-y(1,1:end-1))/(xpixel);
ydmat = repmat(yd',[numberofthetabins,1]);

%%
% xsense = linspace(-1-(sensepixel)/2,1+sensorpixel/2,((length(xin)+1)/subsample));
% xsense = linspace(
% xsense = xrange * transpose(linspace(-1-1/(numberofxbins-1), 1+1/(numberofxbins-1), (0.5*numberofxbins+1)+1));

caustic = zeros(length(xsense),length(z0));

xinmat = repmat(xin,[numberofthetabins,1]);
thetainmat = repelem(thetain, numberofxbins);


%% modeling light field

% plane wave
inten_unit = zeros(numberofthetabins,1);
inten_unit((numberofthetabins+1)/2) = 1;
intensity = repelem(inten_unit,numberofxbins);

% from point source
% intensity = reshape(eye(numberofthetabins),numberofxbins*numberofxbins,1);
% intensity = ones(numberofxbins*numberofxbins,1);

%% refracted by diffuser & propagation to the sensor
thetaout = (n * thetainmat/np + (n/np-1) * ydmat);

xout = repmat(xinmat,[1,length(z0)]) + thetaout * z0;

%% accumulate the beam that falls into the sensor pixel

classify = zeros(size(xout,1), size(xout,2));

for j=1:size(xout,2)
    classify(:,j) = discretize(xout(:,j),sensoredges);
    for i=1:size(xout,1)
        if ~isnan(classify(i,j))
            caustic(classify(i,j),j) = caustic(classify(i,j),j) + intensity(i);
        end
    end
end

caustic = caustic + 0.01*randn(size(caustic,1),size(caustic,2));

%% summarize the result
figure(1)
subplot(2,1,1);
hold on;
plot(Dheight);
plot(y);
legend('dh','y');
subplot(2,1,2);
plot(yd);
legend boxoff;
figure(2)
for i=1:length(z0)
    subplot(length(z0),1,i);
    plot(xsense,caustic(:,i))
    xlim([-1000 1000]);
    ylim([-1 10]);
    title(strcat('F=',string(Fresnel(i))));
end

% xsense = zeros(1,length(xin));

