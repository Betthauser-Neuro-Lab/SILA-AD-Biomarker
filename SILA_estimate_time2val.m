function val = SILA_estimate_time2val(tsila,adtime)
% this function is used to calculate a value at a given time from the
% positivity threshold.

extyrs = 3;
%% Create extrapolated model
md1 = fitlm(tsila.adtime(tsila.adtime>max(tsila.adtime)-extyrs),tsila.val(tsila.adtime>max(tsila.adtime)-extyrs));
md2 = fitlm(tsila.adtime(tsila.adtime<min(tsila.adtime)+extyrs),tsila.val(tsila.adtime<min(tsila.adtime)+extyrs));
slopeu = md1.Coefficients.Estimate(2);
intu = md1.Coefficients.Estimate(1);
slopel = md2.Coefficients.Estimate(2);
intl = md2.Coefficients.Estimate(1);

% resample nonparametric curve to finer grid using 0.01 year spacing
tt = min(tsila.adtime):0.01:max(tsila.adtime);
mval = interp1(tsila.adtime,tsila.val,tt);

val = nan(numel(adtime),1);
for i = 1:numel(adtime)
    % create discretely sampled curve for extrapolated values
    % extrapolate to twice the uppoer and lower modeled values
    if tsila.val(end)>tsila.val(1)
        % for increasing with time
        if adtime(i) < min(tsila.adtime)
            % extraplote backwards
            val(i) = adtime(i)*slopel + intl;            
        elseif adtime(i) > max(tsila.adtime)
            %extrapolate forwards
            val(i) = adtime(i)*slopeu + intu;
        else
            % lookup value
            [~,id] = min(abs(tt - adtime(i)));
            val(i) = mval(id);
        end
    elseif tsila.val(end)<tsila.val(1)
        % for increasing with time
        if adtime(i) < min(tsila.adtime)
            % extraplote backwards
            val(i) = adtime(i)*slopeu + intu;
        elseif adtime(i) > max(tsila.adtime)
            %extrapolate forwards
            val(i) = adtime(i)*slopel + intl;
        else
            % lookup value
            [~,id] = min(abs(tt - adtime(i)));
            val(i) = mval(id);
        end
    end
end
