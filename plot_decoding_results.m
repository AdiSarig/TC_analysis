% plot decoded data
paradigms = {'masking', 'dcf', 'ib', 'test_ib', 'test_dcf', 'test_masking'};

% plot results
for paradigm = paradigms
    num_subs = length(faces_vis_stat.(paradigm{:}).mvpa);
    time_cell = repelem({faces_vis_stat.(paradigm{:}).time}, 1, num_subs);

    mv_plot_result(faces_vis_stat.(paradigm{:}).mvpa, time_cell{:}, 'combine', 'average', 'mask', faces_vis_stat.(paradigm{:}).mask);
%     hold on
%     plot(faces_vis_stat.(paradigm{:}).time, 0.5*faces_vis_stat.(paradigm{:}).mask, '*')
    current_title = ['Face visibility - ' paradigm{:}];
    title(current_title);
    ylim([0.35 1.05]);
    xlabel('time');
    ylabel('auc');
    grid off
    xline(0,'--','LineWidth',1.5);
    axesHandlesToChildObjects = findobj(gca, 'Type', 'line');
    saveas(gcf, fullfile('decoding_figures', [current_title '.tif']));
    
    mv_plot_result(houses_vis_stat.(paradigm{:}).mvpa, time_cell{:}, 'combine', 'average', 'mask', houses_vis_stat.(paradigm{:}).mask);
%     hold on
%     plot(houses_vis_stat.(paradigm{:}).time, 0.5*houses_vis_stat.(paradigm{:}).mask, '*')
    current_title = ['House visibility - ' paradigm{:}];
    title(current_title);
    ylim([0.35 1.05]);
    xlabel('time');
    ylabel('auc');
    grid off
    xline(0,'--','LineWidth',1.5);
    saveas(gcf, fullfile('decoding_figures', [current_title '.tif']));
end