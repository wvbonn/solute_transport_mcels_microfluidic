% W. V. Bonneuil
% KTH Royal Institute of Technology, Stockholm, Sweden
% 10/2023
% ---
% plot live fraction in microfluidic culture.
% this script assumes that the data files contain the live fraction
% ('phi_l'). this requires the script get_transport_measures to have
% been run after the data have been saved from the simulations

clear
close all

CONFINED = [1==0 1==1];
EXPORT = 1==1;

Da_num = [2 3];
n.Da = numel(Da_num);
for i = 1:n.Da
    Da_str{i} = num2str(Da_num(i));
end
mkr = {'o';'+';'x'};
Rd_num = [0.25 0.5 1];
n.Rd = numel(Rd_num);
for i = 1:n.Rd
    Rd_str{i} = num2str(Rd_num(i));
end
Pe_num = [0 1 2 4 10 20 40 100 1000 10000 100000];
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

a = 0.5; % m-cels radius (mm)
c_in = 0.2; % inlet concentration (mol/m^3)
col_live = [1 0.6 0.6]; % color of the "live zone" where S_c is such that necrosis is prevented


%% data loading
for h = 1:numel(CONFINED)
    if ~CONFINED(h)
        fold = 'Data\Unconfined\';
        b = 0; % height of m-cels above bottom wall [mm]
        X_Phi(h,:,:,:) = S_c.*(1-Da_cub/4); % variable against which to plot phi_L
    else
        fold = 'Data\Confined\';
        b = 0.4;
        X_Phi(h,:,:,:) = S_c.*(1-Da_cub/3.4); % variable against which to plot phi_L
    end
   
    for i = 1:n.Da
        for j = 1:n.Rd
            for k = 1:n.Pe
                load([fold 'out_domains\Da_' Da_str{i} '_Rd_' Rd_str{j} '_Pe_' Pe_str{k} '.mat']);
                Phi_L(h,i,j,k) = phi_l;
            end
        end
    end
end


%% plotting

figure('position',[400 50 800 400],'color','w');
tiledlayout(1,2);
labels = get_subplot_labels('a':'z',8);

Da_1 = reshape(Da_cub,[],1);
Pe_1 = reshape(Pe_cub,[],1);

% UNCONFINED
col_1.r = reshape(col.r(1,:,:,:),[],1);
col_1.g = reshape(col.g(1,:,:,:),[],1);
col_1.b = reshape(col.b(1,:,:,:),[],1);
x_Phi = reshape(X_Phi(1,:,:,:),[],1);
phi_L = reshape(Phi_L(1,:,:,:),[],1);
x_cont = logspace(min(log10(x_Phi)),max(log10(x_Phi)),100);
nexttile;
hold on;grid on;
set(gca,'fontsize',14,'fontname','times');
set(gca,'xscale','log')
xlabel('${\rm S_c}(1-{\rm Da}/{\rm Da}_{\infty})$','Interpreter','latex')
xlim([3e-2 1e2])
xticks([1e-1 1 10 1e2])
ylabel('$\Phi_L$','Interpreter','latex')
ylim([0 1])
for i = 1:n.Da
    id = intersect(find(Da_1==Da_num(i)),find(Pe_1~=Pe_num(end-2:end)));
    scatter(x_Phi(id),phi_L(id),40,[col_1.r(id) col_1.g(id) col_1.b(id)],mkr{i},'linew',1.5,'handlevisibility','off');
    scatter(-1,-1,30,0.4*[1 1 1],mkr{i},'linew',1.5,'displayname',['Da = ' Da_str{i}]);
end
for j = 1:n.Rd
    scatter(-1,-1,30,[col.r(1,1,j,1) col.g(1,1,j,1) col.b(1,1,j,1)],'square','filled','displayname',['${\rm R_d} = $ ' Rd_str{j}]);
end
xl = get(gca,'XLim');
fill([5 xl(end) xl(end) 5],[0 0 1 1],col_live,'edgecolor','none','FaceAlpha',0.1,'HandleVisibility','off');
annotation('textbox',[.33 .5 .1 .1],'String','no necrosis','fontsize',13,'interpreter','latex',...
           'EdgeColor','none','backgroundcolor','w','HorizontalAlignment','center');
legend('fontsize',12,'fontname','times','location','southeast','NumColumns',1,'Interpreter','latex');
annotation('textbox',[.06 .86 .1 .1],'string',labels{1},'FontSize',16,'FontName','times','FontWeight','bold',...
    'verticalalignment','bottom','edgecolor','none');


% CONFINED
col_1.r = reshape(col.r(2,:,:,:),[],1);
col_1.g = reshape(col.g(2,:,:,:),[],1);
col_1.b = reshape(col.b(2,:,:,:),[],1);
x_Phi = reshape(X_Phi(2,:,:,:),[],1);
phi_L = reshape(Phi_L(2,:,:,:),[],1);
x_cont = logspace(min(log10(x_Phi)),max(log10(x_Phi)),100);
nexttile;
hold on;grid on;
set(gca,'fontsize',14,'fontname','times');
set(gca,'xscale','log')
xlabel('${\rm S_c}(1-{\rm Da}/{\rm Da}_{\infty})$','Interpreter','latex')
% xlim([3e-2 1e2])
xticks([1e-1 1 10 1e2 1e3 1e4])
ylabel('$\Phi_L$','Interpreter','latex')
ylim([0 1])
for i = 1:n.Da
    id = find(Da_1==Da_num(i));
    scatter(x_Phi(id),phi_L(id),40,[col_1.r(id) col_1.g(id) col_1.b(id)],mkr{i},'linew',1.5,'handlevisibility','off');
    scatter(-1,-1,30,0.4*[1 1 1],mkr{i},'linew',1.5,'displayname',['Da = ' Da_str{i}]);
end
for j = 1:n.Rd
    scatter(-1,-1,30,[col.r(2,1,j,1) col.g(2,1,j,1) col.b(2,1,j,1)],'square','filled','displayname',['${\rm R_d} = $ ' Rd_str{j}]);
end
xl = get(gca,'XLim');
fill([100 xl(end) xl(end) 100],[0 0 1 1],col_live,'edgecolor','none','FaceAlpha',0.1,'HandleVisibility','off');
annotation('textbox',[.785 .55 .1 .1],'String',{'no necrosis';'for Da = 2'},'fontsize',13,'interpreter','latex',...
           'EdgeColor','none','backgroundcolor','w','HorizontalAlignment','center');
legend('fontsize',12,'fontname','times','location','southeast','NumColumns',1,'Interpreter','latex');
annotation('textbox',[.51 .86 .1 .1],'string',labels{2},'FontSize',16,'FontName','times','FontWeight','bold',...
    'verticalalignment','bottom','edgecolor','none');




%% export
if EXPORT
    exportgraphics(gcf,'Figures\live_fraction_flow_lowDa.png','Resolution',300);
end

