function [y, mean_value, max_value] = norm_seq(x)

mean_value = mean(x);
zero_mean_x = x - mean_value;
max_value = max(abs(zero_mean_x));
if max_value > 0
    y = zero_mean_x / max_value;
else
    y = zero_mean_x;
end