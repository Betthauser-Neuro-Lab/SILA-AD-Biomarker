function tout = SILA_estimate_other(tilla,test,age,subid)
% This function can be used to estimate values and times from threshold for
% timepoints other than those used to estimate age of threshold and time to
% threshold. This can be useful for example when wanting to impute the
% value of a biomarker at a time that there was no observed result within a
% reasonable time frame. 

% tout is table of estimates for each subject
% texmod is table with the discrete extrapolated modeled curve

%% Parse the inputs
p = inputParser();
addRequired(p,'tilla');
addRequired(p,'test');
addRequired(p,'age',@(x) isnumeric(x))
addRequired(p,'subid',@(x) or(isnumeric(x),ischar(x)))

parse(p,tilla,test,age,subid)
tilla = p.Results.tilla;
test = p.Results.test;
age = p.Results.age;
subid = p.Results.subid;

%% Quick check to make sure no missing people for estimating values (i.e. need at least one amyloid scan in test)
if nnz(ismember(subid,test.subid)~=1)
    disp({'Input subjects missing data from test table';...
        'The below subjects are missing ILLA input estimates'})
    disp(subid(ismember(subid,test.subid)~=1))
    opt = input('Do you wish to continue with omitting the above subjects? (y/n)','s');
    if strcmpi(opt,'y')
        if any(ismember(subid,test.subid)==1)
            age = age(ismember(subid,test.subid)==1);
            subid = subid(ismember(subid,test.subid)==1);
        else
            disp({'No Remaining Subjects';'Terminating Estimation'})
            return
        end
    else
        disp('Terminating Estimation')
        return
    end
end

%% Create extrapolated model
extyrs = 3;
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
% dtl = ceil((mll - mll/2) / slopel);
% ttl = (0:0.01:dtl)-dtl + min(tt);
% vall = (ttl-min(tt))*slopel  + min(mval);
% 
% dtu = ceil((mul*2 - mul)/slopeu);
% ttu = (0:0.01:dtu) + max(tt);
% valu = (ttu-max(tt))*slopeu + max(mval);
% 
% tt = [ttl(1:end-1),tt,ttu(2:end)];
% mval = [vall(1:end-1),mval,valu(2:end)];

%% Get estimates for value, time to 0, and +/-
tout = table();
tout.subid = subid;
tout.age = age;
tout = sortrows(tout,{'subid','age'});
% tout.estval(:) = nan;
% tout.estextrap(:) = boolean;

% tout.valt0(:) = mval(tt==0);
% tout.ageref(:) = nan;
% tout.dtageref(:) = nan;
% tout.estval(:) = nan;
% tout.estaget0(:) = nan;
% tout.estdtt0(:) = nan;
% tout.estresid(:) = nan;
% three cases, first, last, all
% 1559

for i = 1:height(tout)
%     ids = tout.subid==tout.subid(i);
    tsub = tout(i,:); % get id and age for values to estimate
    tmodsub = test(test.subid==tout.subid(i),:); % get values from ILLA estimates
        
    tout.estdtt0(i) = tmodsub.estdtt0(1) + (tsub.age - tmodsub.age(1)); % duration positive at age for estimate
    [~,id0] = min(abs(tt - tout.estdtt0(i)));
    if id0~=1 && id0<numel(tt) 
        % if modeled timepoint is within the model range then just solve
        % non-parametric function
        tout.estval(i) = mval(id0); % estimated value at age for estimate        
    elseif id0==1
        % if modeled timepoint is lower than modeled range then
        % extrapoloate on bottom end of curve
        
        tout.estval(i) = (tout.estdtt0(i))*slopel + intl;
    else
        % if modeled range is above the upper modeled range then
        % extrapolate on the top end of the curve
        tout.estval(i) = (tout.estdtt0(i))*slopeu + intu;
    end
        
    tout.dtalign(i) = tsub.age - tmodsub.ageref(1); % the time between the estimated age and the reference age for estimate
    % was value in extraploated range?
    tout.estextrap(i) = tout.estdtt0(i) > max(tilla.adtime) | tout.estdtt0(i) < min(tilla.adtime);
    % was value outside of modeled limits?
    tout.estoutrange(i) = tout.estdtt0(i) > max(tt) | tout.estdtt0(i) < min(tt);
    % was the age for prediction within the observed age range for
    % subject's scans? note this only works if all aPET observations are included in the model estimation table
    tout.obsrange(i) = tsub.age >= tmodsub.minage(1) & tsub.age <=tmodsub.maxage(1);
    tout.dtminage(i) = tsub.age - tmodsub.minage(1);
    tout.dtmaxage(i) = tsub.age - tmodsub.maxage(1);
%     
%     
%     switch aevent
%         case 'all'
%             if height(tsub)==1
%                 % put single scan on the curve
%                 [~,id0] = min(abs(mval - tout.val(i)));
%                 tout.ageref(i) = tout.val(i);
%                 tout.dtageref(i) = 0;
%                 tout.estval(i) = mval(id0);
%                 tout.estaget0(i) = tout.age(i) - tt(id0);
%                 tout.estdtt0(i) = tout.age(i) - (tout.age(i) - tt(id0));
%             else
%                 % min(SSQ) for dt to place all scans on curve
%                 idmove = round((tsub.age- tsub.age(1))/0.01); % relative to first scan
%                 ll = idmove+1;ul = numel(tt)-max(idmove);
%                 idts = plus(ll:ul,idmove); % generate all possible query indices
%                 smval = mval(idts); % get modeled values for each query indices
%                 [~,id_opt] = min(sum(minus(smval,tsub.val).^2,1)); % find set that minimizes sum of squared residuals
%                 id_optim = min(idts(:,id_opt)); % get the indices of model for optimized fit
%                 
%                 iddt = round((tout.age(i) - tsub.age(1))/0.01);
%                 idt = id_optim + iddt;
%                 tout.ageref(i) = mean(tsub.age);
%                 tout.dtageref(i) = tout.age(i)-mean(tsub.age);
%                 tout.estval(i) = mval(idt);
%                 tout.estaget0(i) = tsub.age(1) - tt(id_optim);
%                 tout.estdtt0(i) = tout.age(i) - (tsub.age(1) - tt(id_optim));     
%             end
%         case 'first'
%             [~,id0] = min(abs(mval - tsub.val(1)));
%             tout.ageref(i) = tsub.age(1);
%             tout.dtageref(i) = tout.age(i)-tsub.age(1);
%             tout.estval(i) = mval(id0 + round((tout.age(i)-tsub.age(1))/0.01));
%             tout.estaget0(i) = tsub.age(1) - tt(id0);
%             tout.estdtt0(i) = tout.age(i) - (tsub.age(1) - tt(id0));
%         case 'last'
%             [~,id0] = min(abs(mval - tsub.val(end)));
%             tout.ageref(i) = tsub.age(end);
%             tout.dtageref(i) = tout.age(i)-tsub.age(end);
%             tout.estval(i) = mval(id0 + round((tout.age(i)-tsub.age(end))/0.01));
%             tout.estaget0(i) = tsub.age(end) - tt(id0);
%             tout.estdtt0(i) = tout.age(i) - (tsub.age(end) - tt(id0));
%     end
end
% tout.estresid = tout.val - tout.estval;
% tout.estpos = tout.estval>=mval(tt==0);

