% W. V. Bonneuil
% KTH Royal Institute of Technology, Stockholm, Sweden
% 10/2023
% ---
% calculate distance from M-CELS (in quasi-1D domain) where concentration
% is above a certain threshold. The calculation is based on the quasi-1D
% analytical expression of concentration in that domain, with the necrotic
% radius equated to 0

clear
close all
EXPORT = 1==1;

Da = [2 3];
n.Da = numel(Da);
mkr = {'o';'+';'x'};
Rd = [0.25 0.5 1];
n.Rd = numel(Rd);
C = [0.95 0.98 0.99];
n.C = numel(C);

for i = 1:n.Rd
    col.r(i) = 0.5-0.5*(i-1)/(n.Rd-1);
    col.g(i) = 0.95-0.35*(i-1)/(n.Rd-1);
    col.b(i) = 0.5-0.5*(i-1)/(n.Rd-1);
end

% convection needed to prevent necrosis
Sc = 5./(1-Da/4);
Pe = Da'.*Sc'./(2*Rd);
% radius of influence
for i = 1:n.C
    R(i,:,:) = -1./Pe.*log(exp(-Pe) - (C(i)-1)*repmat(Sc',1,3));
end


figure('position',[50 50 400 400],'color','w');
hold on;grid on;
set(gca,'fontname','times','fontsize',14);
xlabel('$\lambda$','Interpreter','latex');
ylabel('$r_{\lambda}/L$','Interpreter','latex');
xlim([0.94 1])
ylim([0 0.25])
xticks(C)
for i = 1:n.Da
   for j = 1:n.Rd
        scatter(C,R(:,i,j),50,[col.r(j) col.g(j) col.b(j)],mkr{i},'LineWidth',1.5,'HandleVisibility','off');
   end
end
for i = 1:n.Da
    scatter(2,1,50,0.4*[1 1 1],mkr{i},'linew',1.5,'displayname',['Da = ' num2str(Da(i))]);
end
for j = 1:n.Rd
    scatter(2,1,50,[col.r(j) col.g(j) col.b(j)],'square','filled','displayname',['${\rm R_d} = $ ' num2str(Rd(j))]);
end
legend('fontsize',12,'fontname','times','location','northwest','NumColumns',1,'Interpreter','latex');

if EXPORT
    exportgraphics(gcf,'Figures\sifig_radius_influence.png','Resolution',300)
end
