function params = gammaParams(eventprobs,lens,shape)
% @ eventprobs: [samples(375)*ncases*nstates] = [375*trials*nBumps]
% @ lens: array of trial lengths
% @ shape: 2
% given that the shape is fixed the calculation of the maximum likelihood
% scales becomes simple.  One just calculates the means expected lengths of
% the flats and divides by the shape=

%width=5; % original
width=6;                        % bump width with sampling frequency of 125 Hz
offset = 3;                     % offset to middle point of the bump, originaly = 2
nstates = size(eventprobs,3);   % number of bumps 

averagepos=cat(2,sum(repmat([1:375]',1,nstates).*reshape(mean(eventprobs,2),375,nstates)),mean(lens));
% concatenate horizontaly
averagepos=averagepos-(offset+[0:width:(nstates-1)*width,(nstates-1)*width+offset]);

flats=averagepos-cat(2,0,averagepos(1:end-1));
params=zeros(nstates+1,2);
params(:,1)=shape;
params(:,2)=flats'/shape;
%correct flats between bumps for the fact that the gamma is calculated at
%midpoint

params(2:end-1,2)=params(2:end-1,2)+ .5/shape;

%first flat is bounded on left while last flat may go beyond on right
% This is a bug! -> params(2:end-1,2)=params(2:end-1,2)- .5/shape;
params(1,2)=params(1,2)-.5/shape;
