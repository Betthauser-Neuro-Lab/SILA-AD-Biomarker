function [tsila,tdrs] = SILA(age,value,subid,dt,val0,maxi,varargin)
% SILLA2 sampled iterative local approximation with smoothing
%   [tsila,tdrs] = SILLA2(age,value,subid,dt,val0,maxi) applies the SILA algorithm
%   to input data age, value, subid to approximate a value vs. time curve.
%   dt is the integration interval for Euler's Method used in the
%   approximation of the integrated curve. val0 specifies the initial
%   condition such that f(t=0) = val0. maxi is the maximum number of
%   iterations allowed before the model terminates integration.
%
%   [tsila,tdrs] = SILLA2(age,value,subid,dt,val0,maxi,sk) is the same as above but
%   allows the user to specifiy the size of the smoothing kernel used to
%   smooth the rate vs. value curve. If this value is unspecified, the
%   algorithm will perform an initial step to select a smoothing kernel 
%   that optimizes backwards prediction residuals based on the first and
%   last observations for each person
% 
%   Input Variables:
%       age = number vector corresponding to the age at observation
%       value = number vector with the observed value to be modeled over time
%       subid = number corresponding to a subject identifier
%       dt = number specifying the step size to use for numerical integration
%       val0 = number for the value that corresponds to t=0
%       maxi = number specifying the maximum number of iterations before
%           the model teriminates
%       varargin = optional input number specifying the size of the
%           smoothing kernel
%   Output Variables:
%       tsila = table with discrete value vs. time curve and additional
%       stats info
%       tdrs = table with the discrete values used as input to the
%       intergation of value vs. time with addtiional stats info
%   
%   See Also ILLA, SILA_estimate.
%       
% Written By: Tobey J Betthauser, PhD
%             Univsersity of Wisconsin-Madison
%             Department of Medicine
%             Division of Geriatrics
%
% This function was initiallty designed to estimate the amyloid vs time curve from
% longitudinal PET imaging data. This function calls a subfunction to
% optimize a smoothing kernel of the sampled rate vs level curve, and then
% outputs the optimized model. 
% If using or further developing this method, please cite the following
% reference: Betthauser, et al,. Multi-method investigation of factors 
% influencing amyloid onset and impairment in three cohorts. Brain, 2022

%% create table object and identifiers for order and number of scans
t = table; % create blank table
t.age = age; % create age variable in table
t.val = value; % create value variable in table
t.subid = subid; % create subject ID variable in table
t = sortrows(t,{'subid','age'}); % sort table by subject ID then age
subs = unique(t.subid); % get all unique subject identifiers
for i = 1:numel(subs)
    ids = t.subid==subs(i);
    ages = t.age(ids);
    [~,ida] = sort(ages);
    t.idx(ids) = ida; % for each subject create ordered observation numbers
    t.ns(ids) = numel(ages); % for each subject get the number of longitudinal observations
end
t = t(t.ns>1,:); % remove cases without longitudinal data

%% Setup variables for A+/- ids and residual weighting
% Emphasis will be on backwards prediction. Residuals weighted such that A+
% and A- have equal say despite possible imbalance in the data.
resnorm = nnz(t.val(t.idx==t.ns)>=val0)/nnz(t.val(t.idx==t.ns)<val0); %ratio of positive:negative cases used to weight residuals when estimating smoothing kernel
idpos = t.val>val0; %indices for biomarker negative cases
idneg = t.val<=val0; %indices for biomarker positive cases

%% Run First iteration to identify the smoothing kernel
switch nargin
    case 7
        % use prespecified smoothing kernel
        % note that the data structure for a loop was maintained, but there
        % is only one iteration of the loop for the case with a
        % pre-specified smoothing kernel
        sk = varargin{1};
        dat = struct();
        SSQpos = nan(numel(sk),1);
        SSQneg = nan(numel(sk),1);
        for i = 1:numel(sk)
            [dat.tilla{i},dat.tdrs{i}] = ILLA(t.age,t.val,t.subid,dt,val0,maxi,sk(i));
            temp = SILA_estimate(dat.tilla{i},t.age,t.val,t.subid,'align_event','last','truncate_aget0','no');
            SSQpos(i) = sum(temp.estresid(idpos).^2);
            SSQneg(i) = sum(temp.estresid(idneg).^2);
        end
        [~,ids] = min(SSQpos + resnorm*SSQneg);
        tsila = dat.tilla{ids};
        tdrs = dat.tdrs{ids};
    case 6
        % optimize smoothing kernel
        sk = 0:0.05:0.5; % these values correspond to setting span of smoothing kernel to 0-50% of the data
        dat = struct();
        SSQpos = nan(numel(sk),1);
        SSQneg = nan(numel(sk),1);
        % for each iteration of the loop, the model is estimated and
        % backwards prediction is performed to optimize residuals
        for i = 1:numel(sk)
            [dat.tilla{i},dat.tdrs{i}] = ILLA(t.age,t.val,t.subid,dt,val0,maxi,sk(i));
            temp = SILA_estimate(dat.tilla{i},t.age,t.val,t.subid,'align_event','last','truncate_aget0','no');
            SSQpos(i) = sum(temp.estresid(idpos).^2);
            SSQneg(i) = sum(temp.estresid(idneg).^2);
        end
        % the model with the best backwards prediction is selected using
        % weighted residuals
        [~,ids] = min(SSQpos + resnorm*SSQneg); % get the index of the best fit
        tsila = dat.tilla{ids}; % output the model table with the best fit
        tdrs = dat.tdrs{ids}; % output the discrete rate table from the best fit
end

% Lines below were used to plot residuals and investigate weighting
% parameters
% subplot(3,1,1),plot(sk,SSQpos + resnorm*SSQneg,'.')
% subplot(3,1,2),plot(sk,SSQpos,'.r')
% subplot(3,1,3),plot(sk,resnorm*SSQneg,'.k')

%% Add a variable to specify the optimal smoothing kernel size to the descrete rate table
tdrs.skern(:) = sk(ids);