% W. V. Bonneuil
% KTH Royal Institute of Technology, Stockholm, Sweden
% 10/2023
% ---
% plot necrotic boundary and necrotic center of maximally-supplied M-CELS

clear
close all

EXPORT = 1==1;

Da_num = [4 8 20 50];
n = numel(Da_num);
for i = 1:n
    Da_str{i} = num2str(Da_num(i));
end

a = 0.5; % m-cels radius (mm)
b = 0.4; % height of m-cels above bottom wall [mm]

%% data loading
fold = 'Data\Confined\';
for i = 1:numel(Da_str)
    load([fold 'out_maxsupply\Da_' Da_str{i} '.mat']);
    x_nb{i} = nec_bnd.x/a;
    y_nb{i} = (nec_bnd.y-b)/a;
    x_nc(i) = nec_centre.x/a;
    y_nc(i) = (nec_centre.y-b)/a;
end

%% plotting

% figure creation
figure('position',[50 50 550 400],'color','w');
tiledlayout(1,7);
nexttile([1 6]);hold on;
set(gca,'fontsize',15,'fontname','times');
set(gca,'visible','off');
xlim([-1.1 1.1]); 
ylim([-b/a 1.15]);
% M-CELS surface
theta = pi/180*linspace(-55,235,200);
x_s = cos(theta);
y_s = sin(theta);
plot(x_s,y_s,'linew',1,'color','k');
text(-0.8,0.75,'$\Gamma$','FontSize',15,'Interpreter','latex');
% wall
line([-1 1],-b/a*[1 1],'color','k');
annotation('textbox',[.15 .08 .1 .1],'string','$\mathcal{W}_0$','fontsize',15,'interpreter','latex','edgecolor','none');
x_w = linspace(0.2,0.75,8);
for i = 1:numel(x_w)
    annotation('line',[x_w(i) x_w(i)-0.02],0.11+[0 -0.04],'color','k');
end
% necrotic boundaries
for i = 1:n
    plot(x_nb{i},y_nb{i},'linew',2.5,'color',0.9*(1-(i/n)^1.5)*[0 1 1]+[1 0 0]);
end
% necrotic centre trajectory
plot(x_nc,y_nc,'linew',1,'color','k');
for i = 1:n
    scatter(x_nc(i),y_nc(i),50,0.9*(1-(i/n)^1.5)*[0 1 1]+[1 0 0],'x','linew',2);
end
% necrotic boundaries caption
annotation('textbox',[.7 .73 .1 .1],'string','$\Gamma_n\,(c_n = 10^{-3}\,c_0)$',...
    'fontsize',15,'Interpreter','latex','FitBoxToText','on','EdgeColor','none','VerticalAlignment','top');
annotation('line',[.82 .93],.72*[1 1]);
for i = 1:n
    annotation('line',[0.79 0.85],(0.2+0.4*i/n)*[1 1]+0.06,'linew',2.5,'color',0.9*(1-(i/n)^1.5)*[0 1 1]+[1 0 0]);
    annotation('textbox',[0.85 0.2+0.4*i/n 0.02 0.1],'string',['${\rm Da} = ' Da_str{i} '$'],...
               'fontsize',13,'Interpreter','latex','FitBoxToText','on','EdgeColor','none','VerticalAlignment','middle');
end


if EXPORT
    exportgraphics(gcf,'Figures\necrotic_boundaries_center_maxsupply.png','Resolution',300);
end
