function [tout,texmod] = SILA_estimate(tsila,age,val,subid,varargin)
% SILA_ESTIMATE(t
%   [tsila,tdrs] = SILA(age,value,subid,dt,val0,maxi) applies the SILA algorithm
%   to input data age, value, subid to approximate a value vs. time curve.
%   dt is the integration interval for Euler's Method used in the
%   approximation of the integrated curve. val0 specifies the initial
%   condition such that f(t=0) = val0. maxi is the maximum number of
%   iterations allowed before the model terminates integration.
%
%   [tsila,tdrs] = SILA(age,value,subid,dt,val0,maxi,sk) is the same as above but
%   allows the user to specifiy the size of the smoothing kernel used to
%   smooth the rate vs. value curve. If this value is unspecified, the
%   algorithm will perform an initial step to select a smoothing kernel 
%   that optimizes backwards prediction residuals based on the first and
%   last observations for each person
% 
%   Input Variables:
%       tsila = table output from SILA.m
%       age = number vector corresponding to the age at observation
%       val = number vector with the observed value to be modeled over time
%       subid = number corresponding to a subject identifier
%       varargin = optional input argument to specify name-value pairs for
%       align_event, extrap_years, and truncate_aget0. 
%
%   Output Variables:
%       tsila = table with discrete value vs. time curve and additional
%       stats info
%       tdrs = table with the discrete values used as input to the
%       intergation of value vs. time with addtiional stats info
%
%   Optional Input Argument name-value pairs
%       'align_event', 'first', 'last', or 'all' specifies which
%       observation(s) within a subject are used as the reference to align
%       subject data to the modeled curve. Default is 'last'
%       'extrap_years', numeric specifies the number of years to use to
%       generate a linear model used to extrapolate beyond the modeled
%       curve. Default value is three years.
%       'truncate_aget0', 'yes' or 'no' specifies whether the model should
%       truncate estimated time from threshold that are lower than the
%       lowest time on the modeled curve. Default is 'yes'
%
% Written By: Tobey J Betthauser, PhD
%             Univsersity of Wisconsin-Madison
%             Department of Medicine
%             Division of Geriatrics
%             tjbetthauser@medicine.wisc.edu

%% Parse the inputs
p = inputParser();
addRequired(p,'tilla');
addRequired(p,'val',@(x) isnumeric(x))
addRequired(p,'age',@(x) isnumeric(x))
addRequired(p,'subid',@(x) or(isnumeric(x),ischar(x)))
addParameter(p,'align_event','last',@(x) contains(x,{'first','last','all'}))
addParameter(p,'extrap_years',3,@(x) isnumeric(x))
addParameter(p,'truncate_aget0','yes',@(x) contains(x,{'y','n'}))

parse(p,tsila,val,age,subid,varargin{:})
tsila = p.Results.tilla;
age = p.Results.age;
val = p.Results.val;
subid = p.Results.subid;
aevent = p.Results.align_event;
extyrs = p.Results.extrap_years;
alit0 = p.Results.truncate_aget0;

%% Create extrapolated model
md1 = fitlm(tsila.adtime(tsila.adtime>max(tsila.adtime)-extyrs),tsila.val(tsila.adtime>max(tsila.adtime)-extyrs));
md2 = fitlm(tsila.adtime(tsila.adtime<min(tsila.adtime)+extyrs),tsila.val(tsila.adtime<min(tsila.adtime)+extyrs));
mll = min(tsila.val);mul = max(tsila.val);
slopeu = md1.Coefficients.Estimate(2);
intu = md1.Coefficients.Estimate(1);
slopel = md2.Coefficients.Estimate(2);
intl = md2.Coefficients.Estimate(1);

% resample nonparametric curve to finer grid using 0.01 year spacing
tt = min(tsila.adtime):0.01:max(tsila.adtime);
mval = interp1(tsila.adtime,tsila.val,tt);

% create discretely sampled curve for extrapolated values
% extrapolate to twice the uppoer and lower modeled values
if strcmp(aevent,'all') && tsila.val(end)>tsila.val(1)
    % for increasing with time
    ttl = min(tt*3):0.01:min(tt);ttl = ttl(1:end-1);
    vall = (ttl-min(tt))*slopel  + min(mval);

    ttu = max(tt):0.01:max(tt)*2;ttu = ttu(2:end);
    valu = (ttu-max(tt))*slopeu + max(mval);

    tt = single([ttl(1:end-1),tt,ttu(2:end)]);
    mval = single([vall(1:end-1),mval,valu(2:end)]);
elseif strcmp(aevent,'all') && tsila.val(end)<tsila.val(1)
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
            if dtt0_maxage<min(tsila.adtime)
                dtshift = min(tsila.adtime) - dtt0_maxage;
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

