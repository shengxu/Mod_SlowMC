% H1 Cross Sections

isoname = 'H1';  % isotope name
T = 300;

% !!!IMPORTANT: energy grid for input xs data are in unit MeV

% scattering xs
xs_scat = load([isoname '_' num2str(T) 'K_2']);
scat_size = size(xs_scat, 1);

% capture xs
xs_capt = load([isoname '_' num2str(T) 'K_102']);
capt_size = size(xs_capt, 1);

% fission xs
xs_fiss = [xs_scat(:,1) zeros(scat_size, 1)];
fiss_size = scat_size;

% fiss_size = capt_size;
% xs_fiss = [xs_capt(:,1) 10*ones(fiss_size, 1)];


% set up free gas thermal scattering variables
A = 1;
sizeN = 10000;
kTfactor = [0.01,0.05,0.1,0.25,0.5,0.75,1,2,5,10,15,20,25,30,40,50,75,100,125,150,200,300];
[thermalcdf,Erat]  = free_gas(A,T,sizeN,kTfactor);

% write out HDF5 file
hdfile = horzcat(isoname,'_',num2str(T),'K.h5');
delete(hdfile);

% capture xs data
h5create(hdfile,'/capt_size',1);
h5write(hdfile,'/capt_size',capt_size);
h5create(hdfile,'/engy_capt',capt_size);
h5write(hdfile,'/engy_capt',xs_capt(:,1));
h5create(hdfile,'/xs_capt',capt_size);
h5write(hdfile,'/xs_capt',xs_capt(:,2));
% scattering xs data
h5create(hdfile,'/scat_size',1);
h5write(hdfile,'/scat_size',scat_size);
h5create(hdfile,'/engy_scat',scat_size);
h5write(hdfile,'/engy_scat',xs_scat(:,1)); %same energy grid xs capt xs
h5create(hdfile,'/xs_scat',scat_size);
h5write(hdfile,'/xs_scat',xs_scat(:,2));
% fission xs data
h5create(hdfile,'/fiss_size',1);
h5write(hdfile,'/fiss_size',fiss_size);
h5create(hdfile,'/engy_fiss',fiss_size);
h5write(hdfile,'/engy_fiss',xs_fiss(:,1));
h5create(hdfile,'/xs_fiss',fiss_size);
h5write(hdfile,'/xs_fiss',xs_fiss(:,2));

% for thermal scattering
h5create(hdfile,'/kTsize',1);
h5write(hdfile,'/kTsize',length(kTfactor));

h5create(hdfile,'/cdfsize',1);
h5write(hdfile,'/cdfsize',length(thermalcdf));

h5create(hdfile,'/kT',length(kTfactor));
h5write(hdfile,'/kT',kTfactor);

h5create(hdfile,'/thermalcdf',length(thermalcdf));
h5write(hdfile,'/thermalcdf',thermalcdf);

h5create(hdfile,'/Erat',[size(Erat,1),size(Erat,2)]);
h5write(hdfile,'/Erat',Erat);

h5create(hdfile,'/cdf_width',1);
h5write(hdfile,'/cdf_width',thermalcdf(2) - thermalcdf(1));
