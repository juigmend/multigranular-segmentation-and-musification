function [pind, pval] = peaksfind(v,thr)
%
% DESCRIPTION:
%      Returns indexes and values of peaks in a vector, over a threshold.
%
% SYNTAX:
%      [pind, pval] = peaksfind(v,thr)
%
% INPUT:
%      v : vector
%      thr : threshold
%
% OUTPUT:
%      pind : peaks' indexes
%      pval : peaks' values in a timeseries of length equal to the input vector's length
%
% VERSION: 15 November 2021
%
% Juan Ignacio Mendoza
% University of Jyv?skyl?

peaks_index_all = zeros(1,length(v));
peaks_index_all(2:end-1) = (diff(sign(diff(v))) == -2);
peaks_index_all = v .* peaks_index_all;
peaks_index_selected = (peaks_index_all > thr);
pind = find( peaks_index_selected );
pval = peaks_index_all .* peaks_index_selected;
pval(pval == 0) = NaN;

end