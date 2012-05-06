% O16 Cross Sections

isoname = 'O16';  % isotope name
T = 300;

% !!!IMPORTANT: energy grid for input xs data are in unit MeV

% scattering xs
xs_scat = load([isoname '_' num2str(T) 'K_2']);
scat_size = size(xs_scat, 1);

% capture xs
xs_capt = [xs_scat(:,1) zeros(scat_size, 1)];
capt_size = scat_size;

% fission xs
xs_fiss = [xs_scat(:,1) zeros(scat_size, 1)];
fiss_size = scat_size;

% fiss_size = capt_size;
% xs_fiss = [xs_capt(:,1) 10*ones(fiss_size, 1)];




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