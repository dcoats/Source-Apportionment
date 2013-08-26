function [post] = posterior(F,lamhat,p_lam)
global y unc; %bring in the global values we're using in the funtion

m = lamhat * F;  %calculate our new predicted values of Y

mu = log(m) - .5 * log( (m .^ 2 + unc .^2)./m.^2 );  %take it to the log scale

varl =  log( (m .^ 2 + unc .^2) ./ m.^2) ;

likeli = (-0.5 * ((log(y) - mu).^2)./varl)  - log(y) - .5.*log(varl);  %compute that lognormal density

like = sum(sum(likeli));  %add it up (using two sums sums the matrix)

post = p_lam + like; %add in the prior density to get the un-normalized posterior density

