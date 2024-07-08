function plot_necrotic_boundaries(tilenb,tilesize,X_i,Y_i,X_nec,Y_nec,Pe,CONFINED,param)

    n.Pe = numel(Pe);
    nexttile(tilenb,tilesize); hold on; 
    set(gca,'fontsize',15,'fontname','times');
    set(gca,'visible','off');
    if CONFINED
        xlim([-1.1 1.1]); 
        text(-0.8,0.75,'$\Gamma$','FontSize',15,'Interpreter','latex');
        line([-1 1],-param.b/param.a*[1 1],'color','k');
        annotation('textbox',[.5 .06 .1 .1],'string','$\mathcal{W}_0$','fontsize',15,'interpreter','latex','edgecolor','none');
        % walls
        x.wall = linspace(0.52,0.77,8);
        for i = 1:numel(x.wall)
            annotation('line',[x.wall(i) x.wall(i)-0.02],0.11+[0 -0.04],'color','k');
        end
    else
        xlim([-1.1 1.1]);
        text(-0.8,0.75,'$\Gamma$','FontSize',15,'Interpreter','latex');
        annotation('line',[.1 .45],.11*[1 1],'color','k','lines','--');
        annotation('textbox',[.1 .06 .1 .1],'string','$\mathcal{S}$','fontsize',15,'interpreter','latex','edgecolor','none');
    end
    ylim([-param.b/param.a 1.15])
    % draw fluid-m-cels interface
    plot(X_i,Y_i,'linew',1,'color','k');
    % plot necrotic boundaries
    for k = 1:n.Pe
        plot(X_nec{k},Y_nec{k},'linew',2.5,'color',0.9*(1-(k/n.Pe)^1.5)*[0 1 1]+[1 0 0]);
    end 
end