% 1/v absorber file

% inputs
xs0 = 2;
E0 = 0.0253*1e-6;

% E range
dE = 0.001;
E = 10.^(log10(1e-11):dE:log10(20.0))';

% get size
sizeE = length(E);

% make cross sections
xs_capt = sqrt(E0./E)*xs0;
xs_scat = zeros(sizeE,1);
xs_fiss = zeros(sizeE,1);



% write out hdf5 file
hdfile = 'v_abs.h5';
delete(hdfile);

% capture xs data
h5create(hdfile,'/capt_size',1);
h5write(hdfile,'/capt_size',capt_size);
h5create(hdfile,'/engy_capt',sizeE);
h5write(hdfile,'/engy_capt',E);
h5create(hdfile,'/xs_capt',sizeE);
h5write(hdfile,'/xs_capt',xs_capt);
% scattering xs data
h5create(hdfile,'/scat_size',1);
h5write(hdfile,'/scat_size',scat_size);
h5create(hdfile,'/engy_scat',sizeE);
h5write(hdfile,'/engy_scat',E); %same energy grid xs capt xs
h5create(hdfile,'/xs_scat',sizeE);
h5write(hdfile,'/xs_scat',xs_scat);
% fission xs data
h5create(hdfile,'/fiss_size',1);
h5write(hdfile,'/fiss_size',fiss_size);
h5create(hdfile,'/engy_fiss',sizeE);
h5write(hdfile,'/engy_fiss',E);
h5create(hdfile,'/xs_fiss',sizeE);
h5write(hdfile,'/xs_fiss',xs_fiss);