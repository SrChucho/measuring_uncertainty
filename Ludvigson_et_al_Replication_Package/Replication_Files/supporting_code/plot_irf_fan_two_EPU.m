function fig = plot_irf_fan_two_EPU(Hi,Lo,lim,ylab)
fig=figure;

shock_min(1) = min(min([
    Hi(1,:,[1,2,3]), Hi(4,:,[1,2,3]),Hi(7,:,[1,2,3])...
    Lo(1,:,[1,2,3]), Lo(4,:,[1,2,3]),Lo(7,:,[1,2,3])]));
shock_min(2) = min(min([
    Hi(2,:,[1,2,3]), Hi(5,:,[1,2,3]),Hi(8,:,[1,2,3])...
    Lo(2,:,[1,2,3]), Lo(5,:,[1,2,3]),Lo(8,:,[1,2,3])]));
shock_min(3) = min(min([
    Hi(3,:,[1,2,3]), Hi(6,:,[1,2,3]),Hi(9,:,[1,2,3])...
    Lo(3,:,[1,2,3]), Lo(6,:,[1,2,3]),Lo(9,:,[1,2,3])]));
shock_max(1) = max(max([
    Hi(1,:,[1,2,3]), Hi(4,:,[1,2,3]),Hi(7,:,[1,2,3])...
    Lo(1,:,[1,2,3]), Lo(4,:,[1,2,3]),Lo(7,:,[1,2,3])]));
shock_max(2) =max( max([
    Hi(2,:,[1,2,3]), Hi(5,:,[1,2,3]),Hi(8,:,[1,2,3])...
    Lo(2,:,[1,2,3]), Lo(5,:,[1,2,3]),Lo(8,:,[1,2,3])]));
shock_max(3) = max(max([
    Hi(3,:,[1,2,3]), Hi(6,:,[1,2,3]),Hi(9,:,[1,2,3])...
    Lo(3,:,[1,2,3]), Lo(6,:,[1,2,3]),Lo(9,:,[1,2,3])]));


for i=1:3;
    
    for j=1:3;
        
        scale = 100;
        
        
        k=(i-1)*3+j;
        fig_sub = subplot(3,3,k);
        
        hold on
        ciplot_h( Lo(k,:,1)*scale,Hi(k,:,1)*scale,[0:lim],'k',0.5);
        hold on
        plot([0:lim],  Lo(k,:,2)*scale,'k--','linewidth',2)
        hold on
        plot([0:lim], Hi(k,:,2)*scale,'k--','linewidth',2)
        hold on
        plot([0:lim], zeros(lim+1,1),'k','Linewidth',1)
        title([ylab{i} ' shock'],'Fontsize',10,'fontweight','bold');
        ylabel([ylab{j}],'Fontsize',10);
        if i >2
            xlabel('Months','Fontsize',9);
        end
        xlim([0,lim])
        
        if j <4
            Lo_1_1 = linspace(shock_min(j)*scale,0,2);
            Lo_1_2 = linspace(0,shock_max(j)*scale,3);
        end
        
        
        set(fig_sub ,'YTick',[Lo_1_1(1:end-1),Lo_1_2])
        set(fig_sub ,'Ylim',[shock_min(j)*scale,shock_max(j)*scale])
              tix=get(gca,'ytick')';
        set(gca,'yticklabel',num2str(tix,'%.1f'))
        
        if i == 2 & j == 2
            leg = legend('EPU Base','EPU News');
            set(leg,'Location','southwest','box','off','Fontsize',9)
        else
            
        end
        
        
        
        
        
    end
end