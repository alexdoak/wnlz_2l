function [] = plotTiled(data,axislabels)

M = numel(data);

t = tiledlayout(M,1);
t.TileSpacing = 'compact';
t.Padding = 'compact';

ax = gobjects(M,1);

black = [0,0,0];
red  = [0.75 0.10 0.10];
blue = [0.00 0.35 0.75];
colours = {black,blue,red}

for p = 1:M

    ax(p) = nexttile;

    for q = 1:numel(data(p).y)
        plot(data(p).x{q},data(p).y{q},'Color',colours{q})
        hold on
        colours{q}
    end

    hold off

    if p < M
        ax(p).XTickLabel = [];
    end

    if axislabels 
    ylabel('$\zeta$','Rotation',0,'FontSize',12)
        if p==M
            xlabel('$x$','FontSize',12)
        end
    end

end

linkaxes(ax, 'x')




end



