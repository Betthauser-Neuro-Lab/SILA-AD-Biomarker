function tout = SILA_estimate_other(tsila,test,age,subid)
% This function can be used to estimate values and times from threshold for
% timepoints other than those used to estimate age of threshold and time to
% threshold. This can be useful for example when wanting to impute the
% value of a biomarker at a time that there was no observed result within a
% reasonable time frame. 

% tout is table of estimates for each subject

%% Parse the inputs
p = inputParser();
addRequired(p,'tilla');
addRequired(p,'test');
addRequired(p,'age',@(x) isnumeric(x))
addRequired(p,'subid',@(x) or(isnumeric(x),ischar(x)))

parse(p,tsila,test,age,subid)
tsila = p.Results.tilla;
test = p.Results.test;
age = p.Results.age;
subid = p.Results.subid;

%% Quick check to make sure no missing people for estimating values (i.e. subject needs to be present in test)
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

%% Get estimates for value, time to 0, and +/-
tout = table();
tout.subid = subid;
tout.age = age;
tout = sortrows(tout,{'subid','age'});


for i = 1:height(tout)
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
    tout.estextrap(i) = tout.estdtt0(i) > max(tsila.adtime) | tout.estdtt0(i) < min(tsila.adtime);
    % was the age for prediction within the observed age range for
    % subject's observations? note this only works if all aPET observations are included in the model estimation table
    tout.obsrange(i) = tsub.age >= tmodsub.minage(1) & tsub.age <=tmodsub.maxage(1);
    tout.dtminage(i) = tsub.age - tmodsub.minage(1);
    tout.dtmaxage(i) = tsub.age - tmodsub.maxage(1);
end

