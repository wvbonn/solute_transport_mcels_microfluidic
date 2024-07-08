% W. V. Bonneuil
% KTH Royal Institute of Technology, Stockholm, Sweden
% 10/2023
% ---
% plot necrotic centre trajectory in fluidic culture.
% this script assumes that the data files contain the necrotic centre
% ('nec_centre'). this requires the script get_transport_measures to have
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
mkr = {'.';'x';'*'};
Rd_num = [0.5 1 2];
n.Rd = numel(Rd_num);
for i = 1:n.Rd
    Rd_str{i} = num2str(Rd_num(i));
end
Pe_num = [0.1 0.5 1 2 4 10 20 40 100 1000];
n.Pe = numel(Pe_num);
for i = 1:n.Pe
    Pe_str{i} = num2str(Pe_num(i));
end

for i = 1:n.Da
    for j = 1:n.Rd
        for k = 1:n.Pe
            S_c(i,j,k) = 2*Rd_num(j)*Pe_num(k)/Da_num(i);
            X_cent(i,j,k) = sqrt(Da_num(i))*S_c(i,j,k); % variable against which to plot the position the necrotic centre
        end
            col.r(1,i,j,1:n.Pe) = 0.5-0.5*(j-1)/(n.Rd-1); % declare plot colours
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
    else
        fold = 'Data\Confined\';
        b = 0.4;
    end
   
    for i = 1:n.Da
        for j = 1:n.Rd
            for k = 1:n.Pe
                load([fold 'out_domains\Da_' Da_str{i} '_Rd_' Rd_str{j} '_Pe_' Pe_str{k} '.mat']);
                x_nec(h,i,j,k) = nec_centre.x;
            end
        end
    end
end


%% plotting

figure('position',[50 50 800 450],'color','w');
tiledlayout(6,2);
labels = get_subplot_labels('a':'z',8);

% UNCONFINED
nexttile([4 1]);hold on;
set(gca,'fontsize',15,'fontname','times');
set(gca,'xscale','log')
xlabel('${\rm Da}^{1/2}\,{\rm S_c}$','Interpreter','latex');
xticks([0.1 1 10 100 1000]);
ylabel('$\overline{x}_n/a$','Interpreter','latex');
ylim([-0.01 0.2])
yticks([0 0.1 0.2])
for i = 1:n.Da
    for j = 1:n.Rd
        plot(reshape(X_cent(i,j,:),[],1),reshape(x_nec(1,i,j,:),[],1)/a,'linew',1.5,...
            'color',[col.r(1,i,j,1) col.g(1,i,j,1) col.b(1,i,j,1)],'marker',mkr{i});
    end
end
annotation('textbox',[.05 .89 .1 .1],'string',labels{1},'FontSize',16,'FontName','times','FontWeight','bold',...
    'verticalalignment','bottom','edgecolor','none');
title('unconfined','interpreter','latex')

% CONFINED
nexttile([4 1]);hold on;
set(gca,'fontsize',15,'fontname','times');
set(gca,'xscale','log')
xlabel('${\rm Da}^{1/2}\,{\rm S_c}$','Interpreter','latex');
xticks([0.1 1 10 100 1000]);
ylabel('$\overline{x}_n/a$','Interpreter','latex');
ylim([-0.01 0.2])
yticks([0 0.1 0.2])
for i = 1:n.Da
    for j = 1:n.Rd
        plot(reshape(X_cent(i,j,:),[],1),reshape(x_nec(2,i,j,:),[],1)/a,'linew',1.5,...
            'color',[col.r(2,i,j,1) col.g(2,i,j,1) col.b(2,i,j,1)],'marker',mkr{i});
    end
end
annotation('textbox',[.51 .89 .1 .1],'string',labels{2},'FontSize',16,'FontName','times','FontWeight','bold',...
    'verticalalignment','bottom','edgecolor','none');
title('confined','interpreter','latex')

% legend, unconfined
nexttile(9,[2 1]);hold on;
xlim([0 1]); ylim([0 1]);
set(gca,'visible','off');
for i = 1:n.Da
    scatter(-1,-1,30,0.4*[1 1 1],mkr{i},'linew',1.5,'displayname',['${\rm Da = }$ ' Da_str{i}]);
end
for j = 1:n.Rd
    scatter(-1,-1,30,[col.r(1,1,j,1) col.g(1,1,j,1) col.b(1,1,j,1)],'square','filled','displayname',['${\rm R_d} = $ ' Rd_str{j}]);
end
legend('fontsize',12,'fontname','times','location','south','NumColumns',2,'Interpreter','latex');
% legend, confined
nexttile(10,[2 1]);hold on;
xlim([0 1]); ylim([0 1]);
set(gca,'visible','off');
for i = 1:n.Da
    scatter(-1,-1,30,0.4*[1 1 1],mkr{i},'linew',1.5,'displayname',['${\rm Da = }$ ' Da_str{i}]);
end
for j = 1:n.Rd
    scatter(-1,-1,30,[col.r(2,1,j,1) col.g(2,1,j,1) col.b(2,1,j,1)],'square','filled','displayname',['${\rm R_d} = $ ' Rd_str{j}]);
end
legend('fontsize',12,'fontname','times','location','south','NumColumns',2,'Interpreter','latex');



%% export
if EXPORT
    exportgraphics(gcf,'Figures\necrotic_center_flow.png','Resolution',300);
end

