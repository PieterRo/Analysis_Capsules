% compare the color SNR to the SNR in m1

SNR_Paolo = m1.SNR;

snrPaolo_V1 = mean(SNR_Paolo(1:512, :), 2, 'omitnan');   % [512 x 1]
bestSNR = bestSNR(:);                                % [512 x 1]
ok = isfinite(bestSNR) & isfinite(snrPaolo_V1);
x = bestSNR(ok);
y = snrPaolo_V1(ok);

figure;
scatter(x, y, 'filled');
xlabel('Best SNR (max over color+window)');
ylabel('SNR\_Paolo (mean over 3 conditions, V1)');
title(sprintf('V1 sites (N = %d)', numel(x)));
grid on;

r = corr(x, y, 'Type','Pearson');
p = polyfit(x, y, 1);
xx = linspace(min(x), max(x), 200);
yy = polyval(p, xx);
hold on;
plot(xx, yy, 'LineWidth', 2);
legend('Sites', sprintf('Fit: y=%.3f x + %.3f, r=%.3f', p(1), p(2), r), 'Location','best');

