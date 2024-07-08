% W. V. Bonneuil
% KTH Royal Institute of Technology, Stockholm, Sweden
% 10/2023
% ---
% declare parameters, run simulations and export their results
% (maximally-supplied M-CELS)


clear
close all
tic

CONFINED = 1==0;
% declare parameters
c_0 = 0.2; % inlet concentration (mol/m^3)
a = 0.5; % M-CELS radius (mm)
Da = 5;
% solute diffusivity in M-CELS
param_in(1,:) = 1e-10*ones(1,numel(Da));
% consumption rate
param_in(2,:) = Da*c_0.*param_in(1,:)/a^2;
% Damköhler
param_in(3,:) = Da;
% save path
if CONFINED
    savepath = 'Data\Confined\';
else
    savepath = 'Data\Unconfined\';
end


for j = 1:size(param_in,2)
    if CONFINED
        MCELS_maxsupply_confined(param_in(:,j));
    else
        MCELS_maxsupply_unconfined(param_in(:,j));
    end
    data = load('out_sim_maxsupply.txt');
    x = data(:,1);
    y = data(:,2);
    c = data(:,3);
    save([savepath 'out_maxsupply\Da_' num2str(param_in(3,j)) '.mat'],'x','y','c');
    
    disp(j)
    toc
end



