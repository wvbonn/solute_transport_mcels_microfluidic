% W. V. Bonneuil
% KTH Royal Institute of Technology, Stockholm, Sweden
% 10/2023
% ---
% plot live fraction in static culture.
% this script assumes that the data files contain the live fraction
% ('phi_l'). this requires the script get_transport_measures to have
% been run after the data have been saved from the simulations

clear
close all

CONFINED = [1==0 1==1];
EXPORT = 1==0;

Da_num = [4 8 20];
n.Da = numel(Da_num);
for i = 1:n.Da
    Da_str{i} = num2str(Da_num(i));
end
mkr = {'o';'+';'x'};
Rd_num = [0.1 0.25 0.5 1 2 4 10 20 100];
n.Rd = numel(Rd_num);
for i = 1:n.Rd
    Rd_str{i} = num2str(Rd_num(i));
end

for i = 1:n.Da
    for j = 1:n.Rd
        Da_sq(i,j) = Da_num(i);
        Rd_sq(i,j) = Rd_num(j);
        S_d(i,j) = 2*Rd_num(j)/Da_num(i); % diffusive supply number
        col.r(1,i,j) = 0.5-0.5*(i-1)/(n.Da-1);
        col.g(1,i,j) = 0.95-0.35*(i-1)/(n.Da-1);
        col.b(1,i,j) = 0.5-0.5*(i-1)/(n.Da-1);
        col.r(2,i,j) = 0.8-0.6*(i-1)/(n.Da-1);
        col.g(2,i,j) = 0.8-0.6*(i-1)/(n.Da-1);
        col.b(2,i,j) = 0.95-0.15*(i-1)/(n.Da-1);
    end
end

a = 0.5; % m-cels radius (mm)
c_in = 0.2; % inlet concentration (mol/m^3)
col_live = [1 0.6 0.6]; % color of the "live zone" where S_d is such that necrosis is prevented


%% data loading
for h = 1:numel(CONFINED)
    if ~CONFINED(h)
        fold = 'Data\Unconfined\';
        b = 0; % height of m-cels above bottom wall [mm]
        X_Phi(h,:,:) = Rd_sq.*S_d; % variable against which to plot phi_L
    else
        fold = 'Data\Confined\';
        b = 0.4;
        X_Phi(h,:,:) = Rd_sq.*S_d; % variable against which to plot phi_L
    end
    
    for i = 1:n.Da
        load([fold 'out_maxsupply\Da_' Da_str{i}],'phi_l');
        Phi_L_inf = phi_l;
        for j = 1:n.Rd
            load([fold 'out_domains\Da_' Da_str{i} '_Rd_' Rd_str{j} '_Pe_0.mat']);  
            Phi_L(h,i,j) = phi_l;
        end
        Phi_L(h,i,:) = Phi_L(h,i,:)/Phi_L_inf;
    end
end


%% plotting

figure('position',[400 50 800 400],'color','w');
tiledlayout(1,2);
labels = get_subplot_labels('a':'z',8);

% UNCONFINED
x_Phi = reshape(X_Phi(1,:,:),[],1);
phi_L = reshape(Phi_L(1,:,:),[],1);
x_cont = logspace(min(log10(x_Phi)),max(log10(x_Phi)),100);
nexttile;
hold on;grid on;
set(gca,'fontsize',14,'fontname','times');
set(gca,'xscale','log')
xlabel('${\rm R_d}\,{\rm S_d}$','Interpreter','latex')
xlim([min(x_Phi) max(x_Phi)])
xticks([1e-2 1 1e2])
ylabel('$\Phi_L^0/\Phi_L^\infty$','Interpreter','latex')
ylim([0 1])
scatter(x_Phi,phi_L,40,...
        [reshape(col.r(1,:,:),[],1) reshape(col.g(1,:,:),[],1) reshape(col.b(1,:,:),[],1)],...
        'x','linew',1.5,'handlevisibility','off');
for i = 1:n.Da
    scatter(-1,-1,30,[col.r(1,i,1) col.g(1,i,1) col.b(1,i,1)],'x','linew',1.5,'displayname',['${\rm Da} = $ ' Da_str{i}]);
end
plot(x_cont,-x_cont.^-0.5+(1+x_cont.^-1).^0.5,'lines','--','color','k','handlevisibility','off');
legend('fontsize',12,'fontname','times','location','southeast','NumColumns',1,'Interpreter','latex');
annotation('textbox',[.06 .86 .1 .1],'string',labels{1},'FontSize',16,'FontName','times','FontWeight','bold',...
    'verticalalignment','bottom','edgecolor','none');


% CONFINED
x_Phi = reshape(X_Phi(2,:,:),[],1);
phi_L = reshape(Phi_L(2,:,:),[],1);
x_cont = logspace(min(log10(x_Phi)),max(log10(x_Phi)),100);
nexttile;
hold on;grid on;
set(gca,'fontsize',14,'fontname','times');
set(gca,'xscale','log')
xlabel('${\rm R_d}\,{\rm S_d}$','Interpreter','latex')
xlim([min(x_Phi) max(x_Phi)])
xticks([1e-2 1 1e2])
ylabel('$\Phi_L^0/\Phi_L^\infty$','Interpreter','latex')
ylim([0 1])
scatter(x_Phi,phi_L,40,...
        [reshape(col.r(2,:,:),[],1) reshape(col.g(2,:,:),[],1) reshape(col.b(2,:,:),[],1)],...
        'x','linew',1.5,'handlevisibility','off');
for i = 1:n.Da
    scatter(-1,-1,30,[col.r(2,i,1) col.g(2,i,1) col.b(2,i,1)],'x','linew',1.5,'displayname',['${\rm Da} = $ ' Da_str{i}]);
end
plot(x_cont,-x_cont.^-0.5+(1+x_cont.^-1).^0.5,'lines','--','color','k','handlevisibility','off');
legend('fontsize',12,'fontname','times','location','southeast','NumColumns',1,'Interpreter','latex');
annotation('textbox',[.51 .86 .1 .1],'string',labels{2},'FontSize',16,'FontName','times','FontWeight','bold',...
    'verticalalignment','bottom','edgecolor','none');



%% export
if EXPORT
    exportgraphics(gcf,'Figures\live_fraction_static_highDa.png','Resolution',300);
end

