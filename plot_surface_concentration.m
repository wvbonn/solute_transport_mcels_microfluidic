function plot_surface_concentration(tilenb,tilesize,theta_i,C_i,Pe_str)   

    n.Pe = numel(Pe_str);
    nexttile(tilenb,tilesize); hold on;
    set(gca,'fontsize',15,'fontname','times');
    xlabel('\theta')
    xlim([0 pi]+pi/180*[-50 50])
    xticks([0 pi/2 pi])
    xticklabels({'\pi';'\pi/2';'0'})
    ylabel('$\gamma$','Interpreter','latex')
    ylim([0 1])
    yticks(0:0.5:1)
    for k = 1:n.Pe
        plot(fliplr(theta_i),C_i{k},'linew',1.5,'color',0.9*(1-(k/n.Pe)^1.5)*[1 1 1]);
    end
    % get the "middle index" of Pe to indicate the corresponding curve on
    % the graph
    id_R = floor(n.Pe/2);
    pos = get(gca,'Position');
    [~,id_C] = min((theta_i-pi).^2); % have the arrow point to \theta=\pi
    annotation('arrow',[pos(1)+0.04 pos(1)+50/280*pos(3)],[pos(2)-0.05 pos(2)+C_i{id_R}(id_C)*pos(4)]);
    annotation('textbox',[pos(1) pos(2)-0.15 .1 .1],'String',['${\rm Pe} = ' Pe_str{id_R} '$'],'FontSize',13,'Interpreter','latex','EdgeColor','none');

end