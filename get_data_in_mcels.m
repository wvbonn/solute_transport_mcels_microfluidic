% W. V. Bonneuil
% KTH Royal Institute of Technology, Stockholm, Sweden
% 10/2023
% ---
% calculate area fraction of M-CELS where solute concentration is above a given threshold

clear
close all

CONFINED = [1==0 1==1];
mkr = {'+';'*';'o'}; % marker types
col = [[0 0.7 0];[0 0 0.8];[0 0 0.8]]; % colors

a = 0.5; % m-cels radius [mm]
b = 0;
fold = 'Data\Unconfined\';

f_unconf = GetFiles([fold 'out_maxsupply\']);
for i = 1:numel(f_unconf)
    load([fold 'out_maxsupply\' f_unconf{i}]); x0 = x; y0 = y; c0 = c;
    X = x0/a;
    Y = (y0-b)/a; % y = 0 at m-cels center, 1 at top, > -1 at bottom 
    id_M = find(X.^2+Y.^2<1);
    x = X(id_M)*a;
    y = Y(id_M)*a+b;
    c = c0(id_M);
    save([fold 'out_maxsupply_M\' f_unconf{i}],'x','y','c');
end

a = 0.5; % m-cels radius [mm]
b = 0.4;
fold = 'Data\Confined\';

f_conf = GetFiles([fold 'out_maxsupply\']);
for i = 1:numel(f_conf)
    load([fold 'out_maxsupply\' f_conf{i}]); x0 = x; y0 = y; c0 = c;
    X = x0/a;
    Y = (y0-b)/a; % y = 0 at m-cels center, 1 at top, > -1 at bottom 
    id_M = find(X.^2+Y.^2<1);
    x = X(id_M)*a;
    y = Y(id_M)*a+b;
    c = c0(id_M);
    save([fold 'out_maxsupply_M\' f_conf{i}],'x','y','c');
end