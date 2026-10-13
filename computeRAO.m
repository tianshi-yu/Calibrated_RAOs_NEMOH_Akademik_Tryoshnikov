%% ComputeRAO calculates the vessel RAO using calibrated hydrodynamic parameters.
%
% mechanics : Calibrated hydrodynamic parameters
% speed     : Vessel speed [knots]
% wavefreq  : Wave angular frequencies [rad/s]
% wavedir   : Wave direction relative to the vessel [deg]

function  RAO = computeRAO(mechanics,speed,wavefreq,wavedir)
    RAO = []; 

    % Indices
    rao_heave = 1; z = 3;
    rao_roll = 2; roll = 4;
    rao_pitch = 3; pitch = 5;

    for dir_i = 1:length(wavedir)
        % dir_i
        Fe = mechanics.Fe{dir_i,1}; % Excitation force
        Inertia = mechanics.Inertia; % Inertia 
        A = mechanics.A{dir_i,1}; % Added mass coeff
        B = mechanics.B{dir_i,1}; % Radiation damping coeff
        Badd = mechanics.Badd{dir_i,1}; % Additional damping 
        Kh = mechanics.Kh; % Stiffness matrix
        Km = mechanics.Km{dir_i,1}; % Additional stiffness matrix
        
        % Ikeda's roll damping
        meshProp = mechanics.meshProp;
        LPP = meshProp.LPP;
        LB = meshProp.LB;
        BD = meshProp.BD;
        CB = meshProp.CB;
        CMID = meshProp.CMID;
        OGD = meshProp.OGD;
        LSKD = meshProp.LSKD;
        BSKB = meshProp.BSKB;
        phi = 20;
        U = speed * 0.514444; %convert knot to m/s
        DispVolume = meshProp.Volume;
        BRollIkeda = calIkedaRollDamping(1./(wavefreq/2/pi),LPP,LB,BD,CB,CMID,OGD,phi,U,DispVolume,LSKD,BSKB,0.1);

        Fe_TF = Fe;
            RAO(:,rao_heave,dir_i) = ( (Fe_TF(:,z)) ./ ...
                (-(Inertia(z,z) + squeeze(A(z,z,:))) .* wavefreq.^2 +...
                1i .* wavefreq .* (squeeze(B(z,z,:))+ Badd(z,z)) + (Kh(z,z) + Km(z,z)) ) );
            RAO(:,rao_roll,dir_i) = ( (Fe_TF(:,roll)) ./ ...
                (-(Inertia(roll,roll) + squeeze(A(roll,roll,:))) .* wavefreq.^2 +...
                1i .* wavefreq .* ( BRollIkeda + Badd(roll,roll)) +...
                (Kh(roll,roll) + Km(roll,roll)) ) );
            RAO(:,rao_pitch,dir_i) = ( (Fe_TF(:,pitch)) ./ ...
                (-(Inertia(pitch,pitch) + squeeze(A(pitch,pitch,:))) .* wavefreq.^2 +...
                1i .* wavefreq .* (squeeze(B(pitch,pitch,:)) + Badd(pitch,pitch)) + (Kh(pitch,pitch) + Km(pitch,pitch)) ) );
            
        
    end
    
end