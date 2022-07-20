function spaghetti_plot(age,val,subid,varargin)
% SPAGHETTI_PLOT generates a spahgetti plot from longitudinal input data.
% Individual subjects with only one observation are plotted as single
% scatter points.
%
% Usage: 
%   spaghetti_plot(age,val,subid): age, val and subid should be 1D numeric
%   vectors with the age at observation and the value to be plotted. subid
%   should also be numeric. All vectors should have the same length.
%
%   spaghetti_plot(age,val,subid,ids): same as above, but will only plot
%   the elements in ids that are set to true. length of ids should be the
%   same as length of age. 
%
% Written By: Tobey J Betthauser, PhD
%             Assistant Professor of Medicine
%             University of Wisconsin-Madison, SMPH
%
% Version History:
%   Created on 20210604

tin = table(age,val,subid,'VariableNames',{'age','val','subid'});
tin = sortrows(tin,{'subid','age'});

if nargin==4
    tin = tin(varargin{1},:);
end

subs = unique(tin.subid);
% figure('Units','inches','Position',[1,1,3,3])
for i = 1:numel(subs)
    ts = tin(tin.subid==subs(i),:);
    plot(ts.age,ts.val,'.-'),hold on
end