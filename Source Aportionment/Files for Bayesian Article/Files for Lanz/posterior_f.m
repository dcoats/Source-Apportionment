function [post] = posterior(F,lamhat,p_F)
global y  unc;  %%%SEE THE COMMENTS ON POSTERIOR.M

m = lamhat * F; %get expectation

mu = log(m) - .5 * log( (m .^ 2 + unc .^2)./m.^2 ); %take to the log scale

st_dev =  log( (m .^ 2 + unc .^2) ./ m.^2) ;

likeli = (-0.5 * ((log(y) - mu).^2)./st_dev)  - log(y) - .5.*log(st_dev); %get density

like = sum(likeli); %here we just sum across the columns instead of the whole thing

post = p_F + like; %add in the prior densityd