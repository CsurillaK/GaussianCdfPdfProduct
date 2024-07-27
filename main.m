% Pre-compute symbolic expression (Requires Symbolic Toolbox)
clearvars

syms m s m1 s1 m2 s2 x t

assume(m, "real");
assume(m1, "real");
assume(m2, "real");
assume(s, "positive");
assume(s1, "positive");
assume(s2, "positive");

normPDF(x) = exp(-x^2/2) / sqrt(2*sym(pi));
normCDF(x) = int(normPDF(t), t, sym(-Inf), x);
gaussPDF(x, m, s) = normPDF((x-m)/s)/s;
gaussCDF(x, m, s) = normCDF((x-m)/s);

p(x) = simplify(gaussPDF(x, m1, s1) * gaussCDF(x, m2, s2) / normCDF((m1-m2)/sqrt(s1^2+s2^2)));

d_log_p(x) = simplify(diff(log(p(x)), x));
dd_log_p(x) = simplify(diff(d_log_p(x), x));
mean_update(m) = simplify(m - d_log_p(m) / dd_log_p(m)); % One step Newton's method
approximate_variance(x) = simplify((-dd_log_p(x))^-0.5);

% Save to MAT-file
p = matlabFunction(p);
mean_update = matlabFunction(mean_update);
approximate_variance = matlabFunction(approximate_variance);
save("./Misc/precomputed_equations.mat", "p", "mean_update", "approximate_variance");

%% Interactive plot
clearvars
addpath("Misc\");
functions = load("./Misc/precomputed_equations.mat");

abscissa = linspace(-1.5, 1.5, 500);
f = figure("Name", "Product of Gaussian PDF and CDF", "NumberTitle", "off", "MenuBar", "none", "ToolBar", "none");
ax = axes();
graphics = struct();
graphics.line_product = line("Parent", ax, "XData", abscissa, "YData", abscissa, "LineWidth", 2, "Color", "k");
graphics.line_pdf = line("Parent", ax, "XData", abscissa, "YData", abscissa, "LineWidth", 1, "Color", "r");
graphics.line_cdf = line("Parent", ax, "XData", abscissa, "YData", abscissa, "LineWidth", 1, "Color", "b");
graphics.line_mean = xline(ax, 0, "--");
graphics.text_mean = text(0, 0, "", "VerticalAlignment", "bottom", "HorizontalAlignment", "left", "FontSize", 10, "Color", "k");
title(ax, "Product of Gaussian PDF and CDF");

manipulate(gcf, @(m1, s1, m2, s2) update_figure(m1, s1, m2, s2, graphics, functions), ...
    {"PDF mean", -1, 1, 0}, {"PDF variance", 0.1, 2, 1}, ...
    {"CDF mean", -1, 1, 0}, {"CDF variance", 0.1, 2, 1}, ...
    "UpdateMode", "High");

legend(ax, [graphics.line_pdf, graphics.line_cdf, graphics.line_product], ["PDF", "CDF", "Product"], "Location", "northwest")
grid(ax);

function update_figure(m1, s1, m2, s2, graphics, functions)
graphics.line_pdf.YData = normpdf(graphics.line_pdf.XData, m1, s1);
graphics.line_cdf.YData = normcdf(graphics.line_cdf.XData, m2, s2);
graphics.line_product.YData = graphics.line_pdf.YData .* graphics.line_cdf.YData / normcdf((m1-m2)/sqrt(s1^2+s2^2));

% Mean update
mean_update_pdf = functions.mean_update(m1, m1, m2, s1, s2);
mean_update_cdf = functions.mean_update(m2, m1, m2, s1, s2);
value_update_pdf = functions.p(mean_update_pdf, m1, m2, s1, s2);
value_update_cdf = functions.p(mean_update_cdf, m1, m2, s1, s2);
if value_update_pdf > value_update_cdf
    line_mean = mean_update_pdf;
    line_mean_value = value_update_pdf;
    line_mean_color = graphics.line_pdf.Color;
else
    line_mean = mean_update_cdf;
    line_mean_value = value_update_cdf;
    line_mean_color = graphics.line_cdf.Color;
end
graphics.line_mean.Value = line_mean;
graphics.line_mean.Color = line_mean_color;
graphics.text_mean.Color = line_mean_color;
graphics.text_mean.Position = [line_mean, line_mean_value];
graphics.text_mean.String = sprintf("Mean: %.4f\nVariance: %.4f", line_mean, functions.approximate_variance(line_mean, m1, m2, s1, s2));
end
