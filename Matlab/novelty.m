function nov = novelty(sm,bw,rs)
%
% DESCRIPTION:
%      Calculates novelty from a distance (A.K.A similarity, dissimilarity) matrix (Foote & Cooper, 2003).
%
% SYNTAX:
%      nov = mcnovelty(sm,bw)
%      nov = mcnovelty(sm,bw,rs)
%
% INPUT:
%      sm : similarity matrix (square matrix)
%      bw : kernel bandwidth (integer)
%      rs : rescale to [1...0], default = 0 (logical)
%
% OUTPUT:
%      nov : a novelty score (ranging from 0 to 1)
%
% VERSION: 15 November 2021
%
% Juan Ignacio Mendoza
% University of Jyv?skyl?

window_width = round( bw * 8);
if rem(window_width,2) % this ensures that the length of the window is even so that the checkerboard can be produced
    window_width = window_width + 1;
end
[kernel_x, kernel_y] = meshgrid((-(window_width-1)/2):((window_width-1)/2), (-(window_width-1)/2):((window_width-1)/2));
variance = bw^2;
height =  1 / ( ( 2 * pi * variance) ) ; % this produces unit volume
kernel = height *  exp( - ( (kernel_x.^2) / ( 2 * variance) ) - ( (kernel_y.^2) / ( 2 * variance) ) );
kernel = kernel .* (kron([-1, 1; 1,  -1],ones(window_width/2)));

amt_w = length(sm);
margin = window_width/2;
beginning = margin + 1;
ending = margin + amt_w;
wis = window_width - 1;
ext = zeros(amt_w + window_width);

% pad beginning and ending corners with averages:
ext(1:window_width,1:window_width) = mean(mean(sm(1:window_width/2,1:window_width/2)));
ext(end-window_width:end,end-window_width:end) = mean(mean(sm(end-window_width/2:end,end-window_width/2:end)));

ext( beginning : ending , beginning : ending ) = sm;
nov = zeros(1, amt_w);

for i_3 = 1:amt_w
    wbeg = i_3;
    wend = wis + i_3 ;
    nov(i_3) = sum(sum(ext( wbeg : wend , wbeg : wend ) .* kernel ));
end

if (nargin > 2) && rs
    minnov = min(nov);
    nov = (nov - minnov) / ( max(nov) - minnov);
end

end
