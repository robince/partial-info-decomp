function Pr=prob(data, nR)

PrPos	= zeros(1,nR-1);
% diff of find of diff trick for counting number of elements
temp 	= sort(data(data>0));	% want vector excluding P(0)
dtemp	= diff([temp;max(temp)+1]);
count 	= diff(find([1;dtemp]));
indx 	= temp(dtemp>0);
PrPos(indx)= count ./ numel(data);	% probability
Pr = [max(1-sum(PrPos),0)  PrPos];
Pr = Pr./sum(Pr);
