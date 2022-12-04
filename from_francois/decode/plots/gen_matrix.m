function im_matrix = gen_matrix(conf, genaroc, type)

nexttile()
im = imagesc([conf.win_disp(1) conf.win_disp(end)],[conf.win_disp(1) conf.win_disp(end)],genaroc); 
axis xy square
hold on;
plot([0 0],ylim(),'w--');
plot(xlim(),[0 0],'w--');
plot(xlim(),ylim(),'w-');
xlabel("Testing time");
ylabel("Training time");
% set(gca,'XTick',0:0.5:1,'XTickLabel',{},'YTick',0:0.5:1,'YTickLabel',{});
% caxis([0.3,0.8]);
colormap('hot');
colorbar;
title([type ' trials'])

im_matrix = im.CData;

end