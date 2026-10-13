clear
%% Calibrated RAO
% RAO{speed,1}(frequency, DoF, direction)
% DoF       : 1 - Heave, 2 - Roll, 3 - Pitch
% speed     : Vessel speed [knots]
% frequency : Wave angular frequencies [rad/s]
% direction : Wave direction relative to the vessel [deg]
load('CalibratedRAO.mat')

%% Compute RAO using calibrated parameters
% Data:
% og_mechanics      : Original parameters obtained from NEMOH
% calib_mechanics   : Calibrated parameters
% speed-list        : Evaluated vessel speeds [knots]
% wavefreq          : Evaluated wave angular frequencies [rad/s]
% wavedir           : Evaluated wave direction relative to the vessel [deg] 
load('CalibratedParameters.mat')

% Compute RAO
RAOCal = [];
for j = 1:length(speed_list)
    RAOCal{j,1} = computeRAO(calib_mechanics,speed_list(j),wavefreq,wavedir);
end
