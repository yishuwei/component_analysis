function [timecourses,topographies] = component_proj(rawdata,projmat,equalize)

% Input ::
% rawdata = Nchannel x Nsample
% projmat = Nchannel x Ncomponent

% Output ::
% timecourses = Ncomponent x Nsample
% topographies =  Ncomponent x Nchannel


if nargin < 3
    equalize = false;
end
if equalize
    C0 = cov(rawdata');
    rawdata = rawdata/sqrt(trace(C0));
end


timecourses = projmat' * rawdata;

rawdata = rawdata - repmat(mean(rawdata,2),1,size(rawdata,2));

topographies = timecourses * rawdata';

end
