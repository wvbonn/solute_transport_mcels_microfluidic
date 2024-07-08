% W. V. Bonneuil
% KTH Royal Institute of Technology, Stockholm, Sweden
% 10/2023
% ---
% plot live fraction at maximal supply.
% this script assumes that the data files contain the live fraction
% ('phi_l'). this requires the script get_transport_measures_maxsupply to
% have been run after the data have been saved from the simulations


clear
close all

CONFINED = [1==0 1==1];
mkr = {'+';'*';'o'}; % marker types
col = [[0 0.7 0];[0 0 0.8];[0 0 0.8]]; % colors
EXPORT = 1==0;

Da_num = {[3:0.2:4.2 5 6];[2:2:20 30:10:60]};
for i = 1:numel(Da_num)
    for j = 1:numel(Da_num{i})
        Da_str{i}{j} = num2str(Da_num{i}(j));
    end
end

a = 0.5; % m-cels radius (mm)
c_in = 0.2; % inlet concentration (mol/m^3)

% coefficients of asymptotic expression phi_L = a*Da^-0.5-b*Da^-1
Da_x = linspace(min(Da_num{2}),max(Da_num{2}),100);
as.a = sqrt(2);
as.b = -1; 
phi_L_as = 2*(as.a.*(Da_x).^-0.5 + as.b*(Da_x).^-1);


%% data loading
for h = 1:numel(CONFINED)
    if ~CONFINED(h)
        fold = 'Data\Unconfined\';
        b = 0;
    else
        fold = 'Data\Confined\';
        b = 0.4; % height of m-cels center above bottom (mm)
    end

    for i = 1:numel(Da_str)
        for ii = 1:numel(Da_str{i})
            % domain
            load([fold 'out_maxsupply\Da_' Da_str{i}{ii} '.mat']);  
            phi_L{i}(h,ii) = phi_l;
        end
    end
end

%% plot f_n

figure('position',[50 50 800 400],'color','w');
tiledlayout(2,2);
labels = get_subplot_labels('a':'z',8);

nexttile(1,[1 1]);
hold on;grid on;
set(gca,'fontsize',14,'fontname','times');
xlabel('Da')
xlim([min(Da_num{1}) max(Da_num{1})])
ylabel('$\Phi_L^\infty$','Interpreter','latex')
ylim([0.75 1])
for h = 1:2
    scatter(Da_num{1},phi_L{1}(h,:),40,col(h,:),mkr{h},'linew',1.5)
end
annotation('textbox',[.06 .87 .1 .1],'string',labels{1},'FontSize',16,'FontName','times','FontWeight','bold',...
    'verticalalignment','bottom','edgecolor','none');
pos(1,:) = get(gca,'Position');

nexttile(2,[2 1]);
hold on;grid on;
set(gca,'fontsize',14,'fontname','times');
xlabel('Da')
xlim([min(Da_num{2}) max(Da_num{2})])
ylabel('$\Phi_L^\infty$','Interpreter','latex')
ylim([0.25 1])
col_fit = [[0.8 0.8 0];[0 0.8 0.8]];
plot(Da_x,phi_L_as,'color',col_fit(1,:),'linew',1.5);
for h = 1:2
    % [f_f,gof_f] = fit(Da_num{2}(2:end)',1-reshape(f_nec{2}(h,2:end),[],1),ft_pow,'StartPoint',[1,-0.5]);
    % plot(Da_x,f_f.a*Da_x.^f_f.b,'color',col_fit(h,:),'linew',1.5);
    scatter(Da_num{2},phi_L{2}(h,:),40,col(h,:),mkr{h},'linew',1.5);
end
annotation('textbox',[.51 .87 .1 .1],'string',labels{2},'FontSize',16,'FontName','times','FontWeight','bold',...
    'verticalalignment','bottom','edgecolor','none');
pos(2,:) = get(gca,'Position');

nexttile(3,[1 1]);
hold on;
set(gca,'visible','off');
xlim([0 1]); ylim([0 1]);
plot([-1 -1],[-1 -1],'color',col_fit(1,:),'linew',1.5);
for h = 1:2
    scatter(-1,-1,40,col(h,:),mkr{h},'linew',1.5);
end
l = legend('unconfined (asymptotic)','unconfined (simulated)','confined (simulated)','fontsize',12,'fontname','times','location','north');
l.Position(2) = l.Position(2)-0.05;


if EXPORT
    exportgraphics(gcf,'Figures\phi_l_maxsupply.png','Resolution',300);
end
