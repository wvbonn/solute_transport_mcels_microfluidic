clear
close all

EXPORT = 1==0;

Da_num = 5;
Da_str = num2str(Da_num);
Rd_num = 0.5;
Rd_str = num2str(Rd_num);
Pe_num = [1 4 10 20 100 10000 1000000];
n.Pe = numel(Pe_num);
for i = 1:n.Pe
    Pe_str{i} = num2str(Pe_num(i));
end
Pe_disp = {'1';'4';'10';'20';'100';'10^4';'10^6'}; % for display in legend

mkr = {'o';'+';'x'};

param.a = 0.5; % m-cels radius (mm)
param.c_in = 0.2; % inlet concentration (mol/m^3)

% figure preparation
figure('position',[50 50 1100 550],'color','w');
t = tiledlayout(9,21);
labels = get_subplot_labels('a':'z',8);

%% Unconfined

fold = 'Data\Unconfined\';
param.b = 0;
for k = 1:n.Pe
    load([fold 'out_interface\Da_' Da_str '_Rd_' Rd_str '_Pe_' Pe_str{k} '.mat']);
    C_i{k} = c/param.c_in;
end
theta_i = get_angular_coordinate(x,y);
X_i = x/param.a;
Y_i = (y-param.b)/param.a;
% surface concentration, unconfined
plot_surface_concentration(1,[3 9],theta_i,C_i,Pe_str);
annotation('textbox',[.08 .89 .1 .1],'string',labels{1},'FontSize',16,'FontName','times','FontWeight','bold',...
    'verticalalignment','bottom','edgecolor','none');

for k = 1:n.Pe
    load([fold 'out_domains\Da_' Da_str '_Rd_' Rd_str '_Pe_' Pe_str{k} '.mat']);
    G_nec_x{k} = nec_bnd.x/param.a;
    G_nec_y{k} = (nec_bnd.y-param.b)/param.a;
    x_nec_c(k) = nec_centre.x/param.a;
    y_nec_c(k) = (nec_centre.y-param.b)/param.a;
end
% necrotic boundaries, unconfined
plot_necrotic_boundaries(5*21+1,[4 9],X_i,Y_i,G_nec_x,G_nec_y,Pe_str,1==0,param);
% necrotic centre trajectory
plot(x_nec_c,y_nec_c,'linew',1,'color','k');
for k = 1:n.Pe
    scatter(x_nec_c(k),y_nec_c(k),50,0.9*(1-(k/n.Pe)^1.5)*[0 1 1]+[1 0 0],'x','linew',2);
end
annotation('textbox',[.08 .35 .1 .1],'string',labels{3},'FontSize',16,'FontName','times','FontWeight','bold',...
    'verticalalignment','bottom','edgecolor','none');


%% Confined

fold = 'Data\Confined\';
param.b = 0.4;

% sort interface points by angle
load([fold 'out_interface\Da_' Da_str '_Rd_' Rd_str '_Pe_0.mat']);
theta_i = get_angular_coordinate(x,y-param.b);
[theta_i,id_sort_i] = sort(theta_i);
X_i = x(id_sort_i)/param.a;
Y_i = (y(id_sort_i)-param.b)/param.a;

for k = 1:n.Pe
    load([fold 'out_interface\Da_' Da_str '_Rd_' Rd_str '_Pe_' Pe_str{k} '.mat']);
    C_i{k} = c(id_sort_i)/param.c_in;
end
% surface concentration, unconfined
plot_surface_concentration(10,[3 9],theta_i,C_i,Pe_str);
annotation('textbox',[.445 .89 .1 .1],'string',labels{2},'FontSize',16,'FontName','times','FontWeight','bold',...
    'verticalalignment','bottom','edgecolor','none');
% surface concentration caption
annotation('arrow',0.8*[1 1],0.92-[0.18 0],'linew',2,'color','k','headstyle','plain');
annotation('textbox',[0.805 1-0.16 0.1 0.02],'String',{'${\rm Pe}$ increasing';['from $' Pe_disp{1} '$ to $' Pe_disp{end} '$']},...
        'fontsize',14,'interpreter','latex','edgecolor','none','FitBoxToText','on');

for k = 1:n.Pe
    load([fold 'out_domains\Da_' Da_str '_Rd_' Rd_str '_Pe_' Pe_str{k} '.mat']);
    G_nec_x{k} = nec_bnd.x/param.a;
    G_nec_y{k} = (nec_bnd.y-param.b)/param.a;
    x_nec_c(k) = nec_centre.x/param.a;
    y_nec_c(k) = (nec_centre.y-param.b)/param.a;
    Phi_L(k) = phi_l;
end
% necrotic boundaries, unconfined
plot_necrotic_boundaries(3*21+10,[6 9],X_i,Y_i,G_nec_x,G_nec_y,Pe_str,1==0,param);
% necrotic centre trajectory
plot(x_nec_c,y_nec_c,'linew',1,'color','k');
for k = 1:n.Pe
    scatter(x_nec_c(k),y_nec_c(k),50,0.9*(1-(k/n.Pe)^1.5)*[0 1 1]+[1 0 0],'x','linew',2);
end
annotation('textbox',[.445 .53 .1 .1],'string',labels{3},'FontSize',16,'FontName','times','FontWeight','bold',...
    'verticalalignment','bottom','edgecolor','none');

% necrotic boundaries caption
annotation('textbox',[.83 .56 .1 .1],'string','$\Gamma_n\,(c_n = K_{1/2})$',...
    'fontsize',15,'Interpreter','latex','FitBoxToText','on','EdgeColor','none','VerticalAlignment','top','HorizontalAlignment','right');
annotation('line',[.82 .92],.58*[1 1]);
for k = 1:n.Pe
    annotation('line',[0.81 0.84],(0.48*k/n.Pe)*[1 1]+0.05,'linew',2.5,'color',0.9*(1-(k/n.Pe)^1.5)*[0 1 1]+[1 0 0]);
    if Phi_L(k)~=1
        annotation('textbox',[0.85 0.48*k/n.Pe 0.02 0.1],'string',['Pe = $' Pe_disp{k} '$'],...
        'fontsize',13,'Interpreter','latex','FitBoxToText','on','EdgeColor','none','VerticalAlignment','middle');
    else
        annotation('textbox',[0.85 0.48*k/n.Pe 0.02 0.1],'string',['Pe = $' Pe_disp{k} '^*$'],...
        'fontsize',13,'Interpreter','latex','FitBoxToText','on','EdgeColor','none','VerticalAlignment','middle');
    end
end

% fluid flow direction arrows
annotation('arrow',[.12 .17],.37*[1 1],'linew',1.5,'headstyle','plain');
annotation('arrow',[.49 .54],.56*[1 1],'linew',1.5,'headstyle','plain');
annotation('textbox',[.11 .35 .1 .1],'string','fluid flow',...
        'fontsize',13,'Interpreter','latex','FitBoxToText','on','EdgeColor','none','VerticalAlignment','middle');
annotation('textbox',[.48 .54 .1 .1],'string','fluid flow',...
        'fontsize',13,'Interpreter','latex','FitBoxToText','on','EdgeColor','none','VerticalAlignment','middle');

if EXPORT
    exportgraphics(gcf,'Figures\necrotic_center_and_boundaries_flow.png','Resolution',600);
end
