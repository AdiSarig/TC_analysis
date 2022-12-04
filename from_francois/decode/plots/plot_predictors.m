function plot_predictors(conf, generalization, time_axis, id)

nt = size(generalization,2);
tic
mean1 = nan(nt,1);
ci1 = nan(nt,2);
for it = 1:nt-1
    mean1(it) = mean(generalization(id.hithigh,it,it));
    ci1(it,:) = bootci(1000, @mean, generalization(id.hithigh,it,it));
end
toc

nexttile(); hold on
plot(time_axis, mean1, 'Color', conf.color.hit)
plot(time_axis, ci1, 'LineStyle', '--', 'Color', conf.color.hit)

end

