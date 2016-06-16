clc
clear all
close all

%% Prepare common parameters
if(~exist('Sat', 'var'))
  Sat = 'C';
end%if

base = '/hpc1/users/Swarm/CalVal4';

%% Ensure paths
if(~exist('Level1bITRF2NEC', 'file'))
  addpath /home/lastec/src/matlab/Swarm/Level1b
end%if
if(~exist('qvxform', 'file'))
  addpath /home/lastec/src/matlab/Swarm/Level1b/quaternions
end%if
if(~exist('xl_sunmoon', 'file'))
  addpath /home/lastec/src/matlab/Swarm/Level1b/EarthExplTools
end%if

%% Prepare variables
t = []; r = []; theta = []; phi = [];
F = []; B = []; 

%% Loop through all files
%files = dir([base '/SW_*_MAG' Sat '_LR_1B*_MDR_MAG_LR.cdf']);
files=dir([base '/SW_*_MAG' Sat '_LR_1B_2015*_MDR_MAG_LR.cdf']);

for filename = { files(:).name }
% filename = {'SW_OPER_MAGC_LR_1B_20150101T000000_20150101T235959_0409_MDR_MAG_LR.cdf'};
% filename = {'SW_OPER_MAGC_LR_1B_20150102T000000_20150102T235959_0409_MDR_MAG_LR.cdf'};

  number = 0;
  % Read required data
%   data = cdfread(fullfile(base, filename{1}), ...
  data = cdfread(fullfile(base, filename{1}), ...
                 'Variables', { 'Timestamp', ...
        'Latitude', 'Longitude', 'Radius', ...
        'F', 'B_VFM', 'q_NEC_CRF', ...
        'Flags_F', 'Flags_B', 'Flags_q', ...
        'dB_Sun' }, ...
                 'ConvertEpochToDatenum', true, 'CombineRecords', true);

  % Select 1-minute, usefull data
  i = find((data{8} < 255) & (data{9} < 255) & (data{10} < 255));
  if Sat=='C'
      i = find((data{9} < 255) & (data{10} < 255));%there's no absolute, so Flag_F
  end
  if(isempty(i))
    continue;   % Just to be safe...
  end%if
  i = i(1:10:end);

  % Append data
  t_tmp = data{1}(i) - datenum(2000,1,1); % MJD2000
  t = [t; t_tmp];
  r = [r; data{4}(i)*1e-3];     % km
  theta = [theta; 90-data{2}(i)];   % deg
  phi = [phi; data{3}(i)];      % deg
  B = [B; data{6}(i,:)];
end

%% Setting F from F_VFM for sat C, for times after november 2014. 
  if Sat=='C' && t(1)>jd2000(2014,11,5)
    F_VFM = sqrt(B(:,1).^2+B(:,2).^2+B(:,3).^2);
    F = [F;F_VFM];
  else
    F = [F; data{5}(i)];
  end
  
%   number = number + 1;
%   percent = number/length(files)
% end%for filename


%% Save .mat file:

%save(sprintf('Sw%c_MagCL_2013_2015.mat',Sat),'t','r','theta','phi','B','F')
save(sprintf('Sw%c_MagCL_2015.mat',Sat),'t','r','theta','phi','B','F')

