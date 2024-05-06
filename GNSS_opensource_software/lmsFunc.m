% Inputs:
%   xn   Input signal sequences
%   dn   Desired signal
%   M    Moving window size
%   mu   step
%   itr  
% Outputs:
%   W    Filter weights  
%   en   Error sequences 
%   yn   Output signal        
function [yn, W]=lmsFunc(xn, dn, N, mu, itr)

en = zeros(itr,1);            
% W  = zeros(N,1);   
W = [0.0113    0.0087    0.0070    0.0051    0.0037    0.0013    0.0061    0.0640    0.9274]';
% Iteration
for epoch = 1:itr                  

    y = W' * xn;        % filter output
    en(epoch) = dn - y ;       % error
    % Weight updating
    W = W + 2*mu*en(epoch)*xn;
end

yn = W'* xn;

end

