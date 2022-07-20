function t = SILA_estimate_val2time(tsila,val)
% SILA_ESTIMATE_VAL2TIME estimates the modeled time from threshold for a
% specified value
%
% t = SILA_ESTIMATE_VAL2TIME(tsila,val) will estimate the time from
% threshold based on the modeled function in tsila for the user specified
% value val. tsila is obtained by from SILA.m
%
% See also SILA, ILLA

%% Parse the inputs
p = inputParser();
addRequired(p,'tilla');
addRequired(p,'val',@(x) isnumeric(x))

parse(p,tsila,val)
tsila = p.Results.tilla;
val = p.Results.val;

%% Create extrapolated model
extyrs = 3;
md1 = fitlm(tsila.adtime(tsila.adtime>max(tsila.adtime)-extyrs),tsila.val(tsila.adtime>max(tsila.adtime)-extyrs));
md2 = fitlm(tsila.adtime(tsila.adtime<min(tsila.adtime)+extyrs),tsila.val(tsila.adtime<min(tsila.adtime)+extyrs));
mll = min(tsila.val);mul = max(tsila.val);
slopeu = md1.Coefficients.Estimate(2);intu = md1.Coefficients.Estimate(1);
slopel = md2.Coefficients.Estimate(2);intl = md2.Coefficients.Estimate(1);

% resample nonparametric curve to finer grid using 0.01 year spacing
tt = min(tsila.adtime):0.01:max(tsila.adtime);
mval = interp1(tsila.adtime,tsila.val,tt);

% create discretely sampled curve for extrapolated values
% extrapolate to twice the uppoer and lower modeled values
% dtl = ceil((mll - mll/2) / slopel);
% ttl = (0:0.01:dtl)-dtl + min(tt);
% vall = (ttl-min(tt))*slopel  + min(mval);
% 
% dtu = ceil((mul*2 - mul)/slopeu);
% ttu = (0:0.01:dtu) + max(tt);
% valu = (ttu-max(tt))*slopeu + max(mval);

% tt = [ttl(1:end-1),tt,ttu(2:end)];
% mval = [vall(1:end-1),mval,valu(2:end)];

%% Get estimates for time to 0, and +/-
slopepos = tsila.val(end)>tsila.val(1);
t = nan(size(val));
for i = 1:numel(val)
    switch slopepos
        case true % increasing value over time
            if val(i)>max(mval)
                t(i) = (val(i) - intu) / slopeu;
            elseif val(i)<min(mval)
                t(i) = (val(i) - intl) / slopel;
            else
                [~,id] = min(abs(val(i)-mval));
                t(i) = tt(id);
            end
        case false % decreasing value over time
            if val(i)<min(mval) 
                t(i) = (val(i) - intu) / slopeu;
            elseif val(i)>max(mval)
                t(i) = (val(i) - intl) / slopel;
            else
                [~,id] = min(abs(val(i)-mval));
                t(i) = tt(id);
            end
    end
end
