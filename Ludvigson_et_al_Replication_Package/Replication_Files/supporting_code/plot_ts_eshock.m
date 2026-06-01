function fig = plot_ts_eshock(e_shock,dates,AA_minus_ind)

%Plot Time Series 
fig = figure;
for k = 1:3
subplot(3,1,k)
shock_plot = zscore(e_shock(k,:))*AA_minus_ind(k);
plot(dates, shock_plot,'linewidth',2);
hold on 
plot(dates, 2*ones(length(dates),1), 'k--', 'Linewidth',1)
hold on 
plot(dates, -2*ones(length(dates),1), 'k--', 'Linewidth',1)
rshade(dates);
xlabel('Year','fontsize',15,'Interpreter','latex')
if k == 1
title('$e_\mathrm{M}$','Fontsize',13,'Interpreter','latex')
elseif k == 2
title('$e_{\mathrm{ip}}$','Fontsize',13,'Interpreter','latex')
else
title('$e_\mathrm{F}$','Fontsize',13,'Interpreter','latex')    
end

% Compute skewness, kurtosis
tx = text(1981,(max(shock_plot))*0.87,sprintf('Skewness = %0.2f and Kurtosis = %0.2f',...
    [skewness(e_shock(k,:));kurtosis(e_shock(k,:))]),'Fontsize',10,'HorizontalAlignment','left','fontweight','bold');
 
%1999 for UFF, 1981 for others
xlim([min(dates),max(dates)])
ylim([min(shock_plot), max(shock_plot)])

end
