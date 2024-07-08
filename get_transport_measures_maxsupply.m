% W. V. Bonneuil
% KTH Royal Institute of Technology, Stockholm, Sweden
% 10/2023
% ---
% get necrotic boundary and its barycentre and calculate area fraction of
% M-CELS where solute concentration is above a given threshold (M-CELS at
% max supply)

clear
close all

CONFINED = [1==0 1==1];
mkr = {'+';'*';'o'}; % marker types
col = [[0 0.7 0];[0 0 0.8];[0 0 0.8]]; % colors

Da_num = {[3:0.2:4.2 5 6];[2:2:40 50:10:100]};
for i = 1:numel(Da_num)
    for j = 1:numel(Da_num{i})
        Da_str{i}{j} = num2str(Da_num{i}(j));
    end
end

a = 0.5; % m-cels radius (mm)
c_in = 0.2; % inlet concentration (mol/m^3)

for h = 1%:numel(CONFINED)
    if ~CONFINED(h)
        fold = 'Data\Unconfined\';
        b = 0;
    else
        fold = 'Data\Confined\';
        b = 0.4; % height of m-cels center above bottom (mm)
    end
    
    % load an interface case to get the m-cels surface area
    load([fold 'out_interface\Da_1_Rd_1_Pe_0.mat']);  
    X_i = x/a;
    Y_i = (y-b)/a;
    theta = get_angular_coordinate(X_i,Y_i);
    [theta,id_sort] = sort(theta);
    A0 = sum((Y_i(id_sort(1:end-1))+b/a).*diff(-X_i(id_sort)));

    for i = 1:numel(Da_str)
        for ii = 1:numel(Da_str{i})
            % load m-cels data
            load([fold 'out_maxsupply\Da_' Da_str{i}{ii} '.mat']);  
            X = x/a;
            Y = (y-b)/a; % y = 0 at m-cels center, 1 at top, > -1 at bottom 
            C_d = c/c_in;
            % contour c = 0.001
            if Da_num{i}(ii)>=10
                id_nec = intersect(find(C_d>0.0005),find(C_d<0.01)); % wider interval at high Da because of high gradients
            else
                id_nec = intersect(find(C_d>0.0005),find(C_d<0.004));
            end            
            % sort points of necrotic boundary by angle
            theta_nec = get_angular_coordinate(X(id_nec),Y(id_nec));
            [theta_nec,id_sortnec] = sort(theta_nec);
            X_nec = X(id_nec(id_sortnec));
            Y_nec = Y(id_nec(id_sortnec));
            % calculate necrotic area
            if ~isempty(id_nec)
                % remove from theta the points where dx is too high
                % (bug in the angle sorting, this solves it)
                id_ang = find(diff(X_nec)>0.1);
                for k = 1:numel(id_ang)
                    X_nec = X_nec([1:id_ang(k)-(k-1) id_ang(k)-(k-1)+2:end]);
                    Y_nec = Y_nec([1:id_ang(k)-(k-1) id_ang(k)-(k-1)+2:end]);
                end
                % if a point is outside the m-cels, replace it by the
                % nearest m-cels-surface point i.e., the one at the same theta
                id_out = find(X_nec.^2+Y_nec.^2>1);
                Y_nec(id_out) = sign(Y_nec(id_out)).*sqrt(abs(1-X_nec(id_out).^2));
                % smooth
                X_nec = smooth(X_nec,5);
                Y_nec = smooth(Y_nec,5);
                % extend necrotic boundary to symmetry/confinement line,
                % which the necrotic points found by intersect don't always
                % meet
                X_nec = [X_nec(1);X_nec;X_nec(end)];
                Y_nec = [-b/a;Y_nec;-b/a];
                A_nec = sum((Y_nec(1:end-1)+b/a).*diff(-X_nec));
                % get coordinates of necrotic core centre: divide necrotic
                % area into small rectangles, weigh the center of those
                % rectangles with their surface area
                da = [Y_nec(1:end-1).*diff(X_nec);(X_nec(end)-X_nec(1))*b/a]; % X_nec may be non-monotonous. The elementary area of the necrotic area is positive if X_nec increases and Y_nec>0 OR if X_nec decreases and Y_nec<0.
                nec_centre.x = sum([X_nec(1:end-1);0]/2.*da)./sum(da)*a;
                if CONFINED(h)
                    nec_centre.y = sum([Y_nec(1:end-1);-b/a]/2.*da)./sum(da)*a+b;
                else
                    nec_centre.y = 0;
                end
                % re-dimensionalise necrotic boundary and pack into structure
                nec_bnd.x = X_nec*a;
                nec_bnd.y = Y_nec*a+b;
            else
                A_nec = 0;
                nec_centre.x = NaN;
                nec_centre.y = NaN;
                nec_bnd.x = [];
                nec_bnd.y = [];
            end
            phi_l = 1-A_nec/A0;
            save([fold 'out_maxsupply\Da_' Da_str{i}{ii} '.mat'],'phi_l','nec_bnd','nec_centre','-append');
            phi_L{i}(h,ii) = 1-A_nec/A0;
        end
    end
end

