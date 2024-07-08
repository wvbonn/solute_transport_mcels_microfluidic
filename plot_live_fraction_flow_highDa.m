% W. V. Bonneuil
% KTH Royal Institute of Technology, Stockholm, Sweden
% 10/2023
% ---
% plot live fraction in fluidic culture.
% this script assumes that the data files contain the live fraction
% ('phi_l'). this requires the script get_transport_measures to have
% been run after the data have been saved from the simulations

clear
close all

CONFINED = [1==0 1==1];
EXPORT = 1==1;

Da_num = [5 10 20];
n.Da = numel(Da_num);
for i = 1:n.Da
    Da_str{i} = num2str(Da_num(i));
end
mkr = {'o';'+';'x'};
Rd_num = [0.5 1 2];
n.Rd = numel(Rd_num);
for i = 1:n.Rd
    Rd_str{i} = num2str(Rd_num(i));
end
Pe_num = [0.1 0.5 1 2 4 10 20 40 100 1000 10000 100000 1000000];
n.Pe = numel(Pe_num);
for i = 1:n.Pe
    Pe_str{i} = num2str(Pe_num(i));
end

for i = 1:n.Da
    for j = 1:n.Rd
        for k = 1:n.Pe
            Da_cub(i,j,k) = Da_num(i);
            Rd_cub(i,j,k) = Rd_num(j);
            Pe_cub(i,j,k) = Pe_num(k);
            S_d(i,j,k) = 2*Rd_num(j)/Da_num(i); % diffusive supply number
            S_c(i,j,k) = 2*Rd_num(j)*Pe_num(k)/Da_num(i); % convective supply number
        end
        col.r(1,i,j,1:n.Pe) = 0.5-0.5*(j-1)/(n.Rd-1);
        col.g(1,i,j,1:n.Pe) = 0.95-0.35*(j-1)/(n.Rd-1);
        col.b(1,i,j,1:n.Pe) = 0.5-0.5*(j-1)/(n.Rd-1);
        col.r(2,i,j,1:n.Pe) = 0.8-0.6*(j-1)/(n.Rd-1);
        col.g(2,i,j,1:n.Pe) = 0.8-0.6*(j-1)/(n.Rd-1);
        col.b(2,i,j,1:n.Pe) = 0.95-0.15*(j-1)/(n.Rd-1);
    end
end
S_ch = S_c.*Da_cub.^0.5; % variable governing concentration asymmetry in the m-cels
tr = [8 20]; % threshold above which concentration is symmetric

a = 0.5; % m-cels radius (mm)
c_in = 0.2; % inlet concentration (mol/m^3)
col_live = [1 0.6 0.6]; % color of the "live zone" where S_d is such that necrosis is prevented


%% data loading
for h = 1:numel(CONFINED)
    if ~CONFINED(h)
        fold = 'Data\Unconfined\';
        b = 0; % height of m-cels above bottom wall [mm]
        X_Phi(h,:,:,:) = Rd_cub.*S_d.*S_c; % variable against which to plot phi_L
    else
        fold = 'Data\Confined\';
        b = 0.4;
        X_Phi(h,:,:,:) = Rd_cub.*S_d.*S_c; % variable against which to plot phi_L
    end
    
    for i = 1:n.Da
        load([fold 'out_maxsupply\Da_' Da_str{i}],'phi_l');
        Phi_L_inf = phi_l;
        for j = 1:n.Rd
            for k = 1:n.Pe
                load([fold 'out_domains\Da_' Da_str{i} '_Rd_' Rd_str{j} '_Pe_' Pe_str{k} '.mat']);
                Phi_L(h,i,j,k) = phi_l/Phi_L_inf;
            end
        end
    end

end


%% plotting

figure('position',[400 50 800 400],'color','w');
tiledlayout(1,2);
labels = get_subplot_labels('a':'z',8);

% UNCONFINED
col_1.r = reshape(col.r(1,:,:,:),n.Da,n.Rd,n.Pe);
col_1.g = reshape(col.g(1,:,:,:),n.Da,n.Rd,n.Pe);
col_1.b = reshape(col.b(1,:,:,:),n.Da,n.Rd,n.Pe);
x_Phi = reshape(X_Phi(1,:,:,:),n.Da,n.Rd,n.Pe);
phi_L = reshape(Phi_L(1,:,:,:),n.Da,n.Rd,n.Pe);
nexttile;
hold on;grid on;
set(gca,'fontsize',14,'fontname','times');
set(gca,'xscale','log')
title(['${\rm Da}^{1/2}\,{\rm S_c}>' num2str(tr(1)) '$'],'Interpreter','latex');
xlabel('${\rm R_d}\,{\rm S_d}\,{\rm S_c}$','Interpreter','latex')
xlim([min(x_Phi(S_ch>tr(1))) max(x_Phi(S_ch>tr(1)))])
xticks([1 1e2 1e4])
ylabel('$\Phi_L/\Phi_L^{\infty}$','Interpreter','latex')
ylim([0 1])
for i = 1:n.Da
    id = intersect(find(S_ch>tr(1)),find(Da_cub==Da_num(i)));
    scatter(x_Phi(id),phi_L(id),40,[col_1.r(id) col_1.g(id) col_1.b(id)],mkr{i},'linew',1.5,'handlevisibility','off');
    scatter(-1,-1,30,0.4*[1 1 1],mkr{i},'linew',1.5,'displayname',['Da = ' Da_str{i}]);
end
for j = 1:n.Rd
    scatter(-1,-1,30,[col.r(1,1,j,1) col.g(1,1,j,1) col.b(1,1,j,1)],'square','filled','displayname',['${\rm R_d} = $ ' Rd_str{j}]);
end
legend('fontsize',12,'fontname','times','location','southeast','NumColumns',1,'Interpreter','latex');
annotation('textbox',[.06 .86 .1 .1],'string',labels{1},'FontSize',16,'FontName','times','FontWeight','bold',...
    'verticalalignment','bottom','edgecolor','none');



% CONFINED
col_1.r = reshape(col.r(2,:,:,:),n.Da,n.Rd,n.Pe);
col_1.g = reshape(col.g(2,:,:,:),n.Da,n.Rd,n.Pe);
col_1.b = reshape(col.b(2,:,:,:),n.Da,n.Rd,n.Pe);
x_Phi = reshape(X_Phi(2,:,:,:),n.Da,n.Rd,n.Pe);
phi_L = reshape(Phi_L(2,:,:,:),n.Da,n.Rd,n.Pe);
nexttile;
hold on;grid on;
set(gca,'fontsize',14,'fontname','times');
set(gca,'xscale','log')
title(['${\rm Da}^{1/2}\,{\rm S_c}>' num2str(tr(2)) '$'],'Interpreter','latex');
xlabel('${\rm R_d}\,{\rm S_d}\,{\rm S_c}$','Interpreter','latex')
xlim([min(x_Phi(S_ch>tr(2))) max(x_Phi(S_ch>tr(2)))])
xticks([1 1e2 1e4])
ylabel('$\Phi_L/\Phi_L^\infty$','Interpreter','latex')
ylim([0 1])
for i = 1:n.Da
    id = intersect(find(S_ch>tr(1)),find(Da_cub==Da_num(i)));
    scatter(x_Phi(id),phi_L(id),40,[col_1.r(id) col_1.g(id) col_1.b(id)],mkr{i},'linew',1.5,'handlevisibility','off');
    scatter(-1,-1,30,0.4*[1 1 1],mkr{i},'linew',1.5,'displayname',['Da = ' Da_str{i}]);
end
for j = 1:n.Rd
    scatter(-1,-1,30,[col.r(1,1,j,1) col.g(1,1,j,1) col.b(1,1,j,1)],'square','filled','displayname',['${\rm R_d} = $ ' Rd_str{j}]);
end
legend('fontsize',12,'fontname','times','location','southeast','NumColumns',1,'Interpreter','latex');
annotation('textbox',[.51 .86 .1 .1],'string',labels{2},'FontSize',16,'FontName','times','FontWeight','bold',...
    'verticalalignment','bottom','edgecolor','none');



%% export
if EXPORT
    exportgraphics(gcf,'Figures\live_fraction_flow_highDa.png','Resolution',300);
end

