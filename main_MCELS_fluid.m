% W. V. Bonneuil
% KTH Royal Institute of Technology, Stockholm, Sweden
% 10/2023
% ---
% declare parameters, run simulations and export their results

clear
close all
tic

CONFINED = 1==0; 
% declare parameters
c_in = 0.2; % inlet concentration [mol/m^3]
Rd = ones(1,2);
Da = 1;
% inlet pressure. Values such that the fluid-channel Péclet number is 1 for
% 3.4e-5 Pa (unconfined) or 2e-5 Pa (confined). Values found empirically
% from simulating fluid flow over the two M-CELS geometries
% param_in(1,:) = repmat([3.4e-6 8.5e-6 1.7e-5 3.4e-5 6.8e-5 1.4e-4 3.4e-4 6.8e-4 1.4e-3 3.4e-3 3.4e-2 3.4e-1 3.4e-1 3.4e-1],1,1);
if CONFINED
    param_in(1,:) = 2e-1*ones(1,numel(Rd));
else
    param_in(1,:) = 3.4e-1*ones(1,numel(Rd));
end
% channel height [mm] (adjusted to keep a gap of 0.3 mm in both confined and
% unconfined configurations)
if CONFINED
    param_in(2,:) = 1.2;
else
    param_in(2,:) = 0.8;
end
% M-CELS radius [mm]
param_in(3,:) = 0.5;
% diffusivity in fluid [m^2/s]
% param_in(4,:) = repmat([1e-9*ones(1,12) 1e-10 1e-11],1,1);
param_in(4,:) = [1e-10 1e-11];
% diffusivity in M-CELS [m^2/s]
param_in(5,:) = param_in(3,:).*param_in(4,:)./(2.5*Rd); % 2.5 = L-a
% consumption [mol/m^3/s]
param_in(6,:) = Da*c_in.*param_in(5,:)./(param_in(3,:)*1e-3).^2;
% Damköhler
param_in(7,:) = Da;
% diffusivity ratio
param_in(8,:) = Rd;
% Péclet
param_in(9,:) = [100000 1000000];
% save path
if CONFINED
    savepath = 'Data\Confined\';
else
    savepath = 'Data\Unconfined\';
end

for j = 1:size(param_in,2)
    if CONFINED
        MCELS_confined(param_in(:,j));
    else
        MCELS_unconfined(param_in(:,j));
    end
    data = load('out_sim_interface.txt');
    x = data(:,1);
    y = data(:,2);
    c = data(:,3);
    save([savepath 'out_interface\Da_' num2str(param_in(7,j)) '_Rd_' num2str(param_in(8,j)) '_Pe_' num2str(param_in(9,j)) '.mat'],'x','y','c');
    
    data = load('out_sim_domains.txt');
    x = data(:,1);
    y = data(:,2);
    c = data(:,3);
    save([savepath 'out_domains\Da_' num2str(param_in(7,j)) '_Rd_' num2str(param_in(8,j)) '_Pe_' num2str(param_in(9,j)) '.mat'],'x','y','c');

    disp(j)
    toc
end

