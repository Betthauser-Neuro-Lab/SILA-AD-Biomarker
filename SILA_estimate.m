function [tout,texmod] = SILA_estimate(tilla,age,val,subid,varargin)

% tout is table of estimates for each subject
% texmod is table with the discrete extrapolated modeled curve

% Remove after testing

%% Parse the inputs
p = inputParser();
addRequired(p,'tilla');
addRequired(p,'val',@(x) isnumeric(x))
addRequired(p,'age',@(x) isnumeric(x))
addRequired(p,'subid',@(x) or(isnumeric(x),ischar(x)))
addParameter(p,'align_event','last',@(x) contains(x,{'first','last','all'}))
addParameter(p,'extrap_years',3,@(x) isnumeric(x))
addParameter(p,'truncate_aget0','yes',@(x) contains(x,{'y','n'}))

parse(p,tilla,val,age,subid,varargin{:})
tilla = p.Results.tilla;
age = p.Results.age;
val = p.Results.val;
subid = p.Results.subid;
aevent = p.Results.align_event;
extyrs = p.Results.extrap_years;
alit0 = p.Results.truncate_aget0;

%% Create extrapolated model
md1 = fitlm(tilla.adtime(tilla.adtime>max(tilla.adtime)-extyrs),tilla.val(tilla.adtime>max(tilla.adtime)-extyrs));
md2 = fitlm(tilla.adtime(tilla.adtime<min(tilla.adtime)+extyrs),tilla.val(tilla.adtime<min(tilla.adtime)+extyrs));
mll = min(tilla.val);mul = max(tilla.val);
slopeu = md1.Coefficients.Estimate(2);
intu = md1.Coefficients.Estimate(1);
slopel = md2.Coefficients.Estimate(2);
intl = md2.Coefficients.Estimate(1);

% resample nonparametric curve to finer grid using 0.01 year spacing
tt = min(tilla.adtime):0.01:max(tilla.adtime);
mval = interp1(tilla.adtime,tilla.val,tt);

% create discretely sampled curve for extrapolated values
% extrapolate to twice the uppoer and lower modeled values
if strcmp(aevent,'all') && tilla.val(end)>tilla.val(1)
    % for increasing with time
    ttl = min(tt*3):0.01:min(tt);ttl = ttl(1:end-1);
    vall = (ttl-min(tt))*slopel  + min(mval);

    ttu = max(tt):0.01:max(tt)*2;ttu = ttu(2:end);
    valu = (ttu-max(tt))*slopeu + max(mval);

    tt = single([ttl(1:end-1),tt,ttu(2:end)]);
    mval = single([vall(1:end-1),mval,valu(2:end)]);
elseif strcmp(aevent,'all') && tilla.val(end)<tilla.val(1)
    % for decreasing with time
    ttl = min(tt*3):0.01:min(tt);ttl = ttl(1:end-1);
    vall = (ttl-min(tt))*slopeu  + max(mval);

    ttu = max(tt):0.01:max(tt)*2;ttu = ttu(2:end);
    valu = (ttu-max(tt))*slopel + min(mval);

    tt = single([ttl(1:end-1),tt,ttu(2:end)]);
    mval = single([vall(1:end-1),mval,valu(2:end)]);
end

%% Get estimates for value, time to 0, and +/-
tout = table();
tout.subid = subid;
tout.age = age;
tout.val = val;
tout = sortrows(tout,{'subid','age'});
tout.minage(:) = nan;
tout.maxage(:) = nan;
tout.valt0(:) = mval(tt==0);
tout.ageref(:) = nan;
tout.dtageref(:) = nan;
tout.estval(:) = nan;
tout.estaget0(:) = nan;
tout.estdtt0(:) = nan;
tout.estresid(:) = nan;

% three cases, first, last, all
for i = 1:height(tout)
    ids = tout.subid==tout.subid(i);
    tsub = tout(ids,:);
    tout.minage(i) = min(tsub.age);
    tout.maxage(i) = max(tsub.age);
    switch aevent
        case 'all'
            if height(tsub)==1
                % put single scan on the curve
                [~,id0] = min(abs(mval - tout.val(i)));
                tout.ageref(i) = tsub.age(1);
                tout.dtageref(i) = 0;
                tout.estval(i) = mval(id0);
                tout.estaget0(i) = tout.age(i) - tt(id0);
                tout.estdtt0(i) = tout.age(i) - (tout.age(i) - tt(id0));
            else
                % min(SSQ) for dt to place all scans on curve
                idmove = round((tsub.age- tsub.age(1))/0.01); % relative to first scan
                ll = idmove+1;ul = numel(tt)-max(idmove);
                idts = plus(ll:ul,idmove); % generate all possible query indices
                smval = mval(idts); % get modeled values for each query indices
                [~,id_opt] = min(sum(minus(smval,tsub.val).^2,1)); % find set that minimizes sum of squared residuals
                id_optim = min(idts(:,id_opt)); % get the indices of model for optimized fit
                
                iddt = round((tout.age(i) - tsub.age(1))/0.01);
                idt = id_optim + iddt;
                tout.ageref(i) = mean(tsub.age);
                tout.dtageref(i) = tout.age(i)-mean(tsub.age);
                tout.estval(i) = mval(idt);
                tout.estaget0(i) = tsub.age(1) - tt(id_optim);
                tout.estdtt0(i) = tout.age(i) - (tsub.age(1) - tt(id_optim));     
            end
        case 'first'
            [~,id0] = min(abs(mval - tsub.val(1)));
            tout.ageref(i) = tsub.age(1);
            tout.dtageref(i) = tout.age(i)-tsub.age(1);
            if id0>=numel(mval) || id0 + round((tout.age(i)-tsub.age(1))/0.01)>numel(mval)
                %extrapolate at the top
                chronfirst = (tsub.val(1)- intu)/slopeu;
                tout.estval(i) = (chronfirst + tout.dtageref(i))*slopeu + intu;
                tout.estdtt0(i) = chronfirst + tout.dtageref(i);
                tout.estaget0(i) = tsub.age(1) - chronfirst;
            elseif id0==1
                % extrapolate from the bottom
                chronfirst = (tsub.val(1)- intl)/slopel;
                tout.estval(i) = (chronfirst + tout.dtageref(i))*slopel + intl;
                tout.estdtt0(i) = chronfirst + tout.dtageref(i);
                tout.estaget0(i) = tsub.age(1) - chronfirst;
            else
                tout.estval(i) = mval(id0 + round((tout.age(i)-tsub.age(1))/0.01));
                tout.estaget0(i) = tsub.age(1) - tt(id0);
                tout.estdtt0(i) = tout.age(i) - (tsub.age(1) - tt(id0));
            end
        case 'last'
            [~,id0] = min(abs(mval - tsub.val(end)));
            tout.ageref(i) = tsub.age(end);
            tout.dtageref(i) = tout.age(i)-tsub.age(end);
            if id0 + round((tout.age(i)-tsub.age(end))/0.01)<1 || id0==1 % this was added 02/05/2021 to fix problems with extrapolation
                chronend = (tsub.val(end)- intl)/slopel;
                tout.estval(i) = (chronend + tout.dtageref(i))*slopel + intl;
                tout.estdtt0(i) = chronend + tout.dtageref(i);
                tout.estaget0(i) = tsub.age(end) - chronend;
            elseif id0>=numel(mval)
                chronend = (tsub.val(end)- intu)/slopeu;
                tout.estval(i) = (chronend + tout.dtageref(i))*slopeu + intu;
                tout.estdtt0(i) = chronend + tout.dtageref(i);
                tout.estaget0(i) = tsub.age(end) - chronend;
            else
                tout.estval(i) = mval(id0 + round((tout.age(i)-tsub.age(end))/0.01));
                tout.estaget0(i) = tsub.age(end) - tt(id0);
                tout.estdtt0(i) = tout.age(i) - (tsub.age(end) - tt(id0));
            end
    end
end
tout.estresid = tout.val - tout.estval;
tout.estpos = tout.estval>=mval(tt==0);
tout.aevent(:) = {aevent};
tout.extrapyrs(:) = extyrs;

%% This is an option to restrict the estaget0 such that at a person's oldest
%   age they cannot have a estaget0 less than minimum ILLA value time to A+

switch contains(alit0,'y')
    case true
        for i = 1:height(tout)
            tsub = tout(tout.subid==tout.subid(i),:);
            dtt0_maxage = tsub.estdtt0(tsub.age==tsub.maxage);
            if dtt0_maxage<min(tilla.adtime)
                dtshift = min(tilla.adtime) - dtt0_maxage;
                tout.estaget0(i) = tout.estaget0(i) - dtshift;
                tout.estdtt0(i) = tout.estdtt0(i) + dtshift;
                tout.truncated(i) = 1;
            else
                tout.truncated(i) = 0;
            end
        end
    case false
        tout.realigned(i) = 0;
end

