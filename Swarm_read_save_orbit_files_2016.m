%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%                                           %%%%%
%%%%% Creating orbit files for Swarm data                 %%%%%
%%%%% Using SwXXX_MagCL_XXXXX_XXXXX.mat as input            %%%%%
%%%%%                                           %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all

addpath /home/nio/m/tools;

%Setting satellite:
Sat ='C';


%% Loading datafile
%file_dir = '/r4/ftp/data/magnetic-satellites/Swarm/SCARF/ASM_VFM_Task_Force/TDS-3/';
%filename_in = sprintf('Sw%c_MagCL_2013_2015.mat',Sat); 
filename_in = sprintf('Sw%c_MagCL_2015.mat',Sat); 

%load([file_dir filename_in])  
load(filename_in)
%t = t - datenum(2000, 1, 1); %correcting time stamp, since the time is given in Matlab time reference, and not MJD2000 as needed

% loads CHAOS-5 model
load '/home/nio/m/mmod/CHAOS-5.mat'

% External field correction (Dst):
filename_Dst = '/r9/cfinl/m/fit/RC_1997-2016_Feb_v1.dat';

[t_Dst, Dst_all, Dst_e_all, Dst_i_all]  = textread(filename_Dst, '%f %f %f %f %*s', 'commentstyle', 'shell');

%% model predictions

% % Selecting only part of file:
% N = 602000; % There are 86000 observations per day. 
% t = t(1:N); r = r(1:N); 
% theta = theta(1:N); 
% phi = phi(1:N); 
% F = F(1:N); 

B_core = synth_values(r, theta, phi, pp, t);
% high degree field
B_crust = synth_values(r, theta, phi, g(1:60*62));
B_int_mod = B_core + B_crust;
% external field (Dst)
RC_ei = interp1(t_Dst, [Dst_e_all Dst_i_all], t, 'linear'); %interpolates the external corrections to the same times as the data is taken.
%B_ext_mod  = synth_values_ext(t, r, theta, phi, m_sm(1:4), m_gsm(1), m_Dst(1), [Dst_e_all Dst_i_all] );
B_ext_mod  = synth_values_ext(t, r, theta, phi, model_ext.m_sm(1:4), model_ext.m_gsm(1), model_ext.m_Dst(1), RC_ei);

a = 6371.2; %radius og the Earth
rad = pi/180;
[QD_lat, QD_lon, Apex_lat, MLT] = qdipole(t, r/a, theta*rad, phi*rad, '/home/nio/m/tools/apexsh_1980-2020.txt');

% loading file with timestamps for orbit shift
if Sat=='A';
tmp = load('/hpc2/users/Swarm/IOV/Auxiliary/Sat_A/AUXAORBCNT/SW_OPER_AUXAORBCNT_20131122T000000_20160204T000000_0001.TXT');
elseif Sat=='B'
tmp = load('/hpc2/users/Swarm/IOV/Auxiliary/Sat_B/AUXBORBCNT/SW_OPER_AUXBORBCNT_20131122T000000_20160204T000000_0001.TXT');
elseif Sat=='C'
tmp = load('/hpc2/users/Swarm/IOV/Auxiliary/Sat_C/AUXCORBCNT/SW_OPER_AUXCORBCNT_20131122T000000_20160204T000000_0001.TXT');
end 

t_orbit = tmp(:,2); % time in MJ2000
n_orbit = tmp(:,1); % gives the orbit number
phi_AN = tmp(:,4); % gives the phi coordinate of the equator crossing

% %Udregn QD_lat og QD_lon: kan evt findes sammen med Apex lat og MLT ved at bruge qdipole(), som her: 
% MLT(1:length(t))=NaN;
% QD_lat(1:length(t))=NaN;
% QD_lon(1:length(t))=NaN;
% Apex_lat=Phi_AN;

%% Save data to orbit files

%lav loop

%lav if loop

%for i:length(t)
%print header
%while t<t_orbit(i+1)
%print værdier: t, r, co-theta, phi, QD_lat, QD_lon, Apex lat, MLT, F, CHAOS(Br,Btheta,Bphi), Ext(Br,Btheta,Bphi)
%end
%t_orbit=t_orbit+1
%end

i_index=1;
i_index_old=0;
filename_out = 'dummy';
for i=1:length(t);
%for i=1:100;
  while t(i) > t_orbit(i_index)
     i_index = i_index + 1;
  end
  if i_index > i_index_old
     fclose('all');
%     [s,err] = unix(['gzip -f ' filename_out]);
%     if s > 0; disp(['Error gzipping file ' filename_out]); end;
%     filename_out = ['/home/cda/Documents/MATLAB/Swarm_orbit_files/',sprintf('%c/%5.5d.abs',Sat,n_orbit(i_index-1))]; % gives the destination and name of the orbit files
     filename_out = [sprintf('/home/jpbrsb/Documents/Orbit%cFiles/',Sat),sprintf('%c_abs/%5.5d.abs',Sat,n_orbit(i_index-1))]; % filename_out = [sprintf('/home/jpbrsb/Documents/Orbit%cFiles/%c_abs',Sat,Sat),sprintf('/%5.5d.abs',Sat,n_orbit(i_index-1))] % gives the destination and name of the orbit files%fullfile('/home/jpbrsb/Documents/OrbitAFiles/',sprintf('%c/%5.5d.vec',Sat,n_orbit(i_index-1))); % gives the destination and name of the orbit files
     fid_out = fopen(filename_out,'w');
     fprintf(fid_out, '%s%s\n', '% SWARM absolute data, CHAOS field model, file created on: ', datestr(now)); % SWARM absolute data
     fprintf(fid_out, '%s\n', '%    t(MJD2000)      r        theta     longitude    QD_lat    QD_lon   Apex_lat   MLT       F        B_r    B_theta     B_phi        B_r    B_theta     B_phi');
     fprintf(fid_out, '%s\n', '%                   [km]      [deg]        [deg]      [deg]     [deg]    [deg]    [hrs]    [nT]   |         CHAOS-model         |          magnetosph.         |');
     fprintf(fid_out, '%s\n', '% ');
     fprintf(fid_out, '%s\n', '% ');
     i_index_old = i_index;
  end
  fprintf(fid_out,'%15.8f %9.3f %11.5f %11.5f %9.3f %9.3f %9.3f %7.3f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f\n',...
  t(i), r(i), theta(i), phi(i), QD_lat(i), QD_lon(i), Apex_lat(i), MLT(i), ...
  F(i), B_int_mod(i,:), B_ext_mod(i,:));
%fclose(fid_out);
end
fclose(fid_out);

%%
% filename_in = '/home/cda/Documents/MATLAB/Swarm_orbit_files/00124.vec';
% tmp = load(filename_in);
% F_obs = tmp(:,9);
% QD_lat = tmp(:,5);
% F_mod = sqrt(tmp(:,10).^2 + tmp(:,11).^2 + tmp(:,12).^2);
% 
% plot(QD_lat, F_obs-F_mod, '.')
