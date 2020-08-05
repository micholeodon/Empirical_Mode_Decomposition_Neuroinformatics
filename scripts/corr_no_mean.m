function c = corr_no_mean(x, y)

norm_x = nansum(abs(x).^2)^(1/2);
norm_y = nansum(abs(y).^2)^(1/2);
c = nansum(x.*y)/(norm_x*norm_y);
