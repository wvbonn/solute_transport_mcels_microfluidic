% W. V. Bonneuil
% KTH Royal Institute of Technology, Stockholm, Sweden
% 10/2023
% ---
% get necrotic boundary and its barycentre and calculate area fraction of
% M-CELS where solute concentration is above a given threshold (M-CELS in
% culture)

clear
close all

CONFINED = [1==0 1==1];

a = 0.5; % m-cels radius [mm]
c_in = 0.2; % inlet solute concentration [mol/m^3]
for h = 1:numel(CONFINED)
    if ~CONFINED(h)
        fold = 'Data\Unconfined\';
        b = 0; % height of m-cels above bottom wall [mm]
    else
        fold = 'Data\Confined\';
        b = 0.4;
    end
    
    % get m-cels area
    load([fold 'out_interface\Da_1_Rd_1_Pe_0.mat']);  
    X_i = x/a;
    Y_i = (y-b)/a;
    theta = get_angular_coordinate(X_i,Y_i);
    [theta,id_sort] = sort(theta);
    A0 = sum((Y_i(id_sort(1:end-1))+b/a).*diff(-X_i(id_sort)));

    files = GetFiles([fold 'out_domains\']);
    for i = 1:numel(files)
        load([fold 'out_domains\' files{i}]);
        k = strfind(files{i},'_');
        p = strfind(files{i},'.');
        Da = str2double(files{i}(k(1)+1:k(2)-1));
        Rd = str2double(files{i}(k(3)+1:k(4)-1));
        Pe = str2double(files{i}(k(5)+1:p(end)-1));
        X = x/a;
        Y = (y-b)/a;
        C_d = c/c_in;
        % contour c = 0.001
        if Rd<0.5
                id_nec = intersect(find(C_d>0.0005),find(C_d<0.002));
        else if Da<10
                id_nec = intersect(find(C_d>0.0005),find(C_d<0.004));
             else
                id_nec = intersect(find(C_d>0.0005),find(C_d<0.01));
             end
        end          
        % sort points of necrotic boundary by angle
        theta_nec = get_angular_coordinate(x(id_nec),y(id_nec));
        [theta_nec,id_sortnec] = sort(theta_nec,'descend');
        X_nec = X(id_nec(id_sortnec));
        Y_nec = Y(id_nec(id_sortnec));
        % calculate necrotic area
        if ~isempty(id_nec)
            % remove from theta the points where dx is too high
            % (bug in the angle sorting, this solves it)
            id_ang = find(abs(diff(X_nec))>0.1);
            for k = 1:numel(id_ang)
                X_nec = X_nec([1:id_ang(k)-(k-1) id_ang(k)-(k-1)+2:end]);
                Y_nec = Y_nec([1:id_ang(k)-(k-1) id_ang(k)-(k-1)+2:end]);
            end
            % if a point is outside the m-cels, replace it by the
            % nearest m-cels-surface point i.e., the one at the same theta
            id_out = find(X_nec.^2+Y_nec.^2>1);
            for l = 1:numel(id_out)
                theta_l = get_angular_coordinate(X_nec(id_out(l)),Y_nec(id_out(l)));
                X_nec(id_out(l)) = cos(theta_l);
                Y_nec(id_out(l)) = sin(theta_l);
            end
            % Y_nec(id_out) = sign(Y_nec(id_out)).*sqrt(abs(1-X_nec(id_out).^2));
            % smooth
            X_nec = smooth(X_nec,5);
            Y_nec = smooth(Y_nec,5);
            % extend necrotic boundary to symmetry/confinement line,
            % which the necrotic points found by intersect don't always
            % meet
            X_nec = [X_nec(1);X_nec;X_nec(end)];
            Y_nec = [-b/a;Y_nec;-b/a];
            A_nec = sum((Y_nec(1:end-1)+b/a).*diff(X_nec));
            % get coordinates of necrotic core centre: divide necrotic
            % area into small rectangles, weigh the center of those
            % rectangles with their surface area
            da = [Y_nec(1:end-1).*diff(X_nec);(X_nec(end)-X_nec(1))*b/a]; % X_nec may be non-monotonous. The elementary area of the necrotic area is positive if X_nec increases and Y_nec>0 OR if X_nec decreases and Y_nec<0.
            nec_centre.x = sum([X_nec(1:end-1);0].*da)./sum(da)*a;
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
        save([fold 'out_domains\' files{i}],'phi_l','nec_bnd','nec_centre','-append');
    end
end

