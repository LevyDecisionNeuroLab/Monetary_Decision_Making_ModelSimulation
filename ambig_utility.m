function y = ambig_utility(base,v,p,AL,alpha,beta,model);

if (strcmp(model,'ambiguity') || strcmp(model,'ambigNrisk')) || strcmp(model,'ambigNriskFixSlope')
    % the model we are using
    y = (p - beta .* (AL./2)) .* v .^alpha + (1-p - beta .* (AL./2)) .* base .^alpha;
elseif strcmp(model,'ambigPower')
    y = p .^ (1+beta.*AL) .* v .^alpha; % change that
elseif strcmp(model,'discounting')
    %y = v ./ (1 + alpha.*log(1+(1-p+beta.*AL./2)./(p-beta.*AL./2)));
    y = v ./ (1 + alpha.*(1-p+beta.*AL./2)./(p-beta.*AL./2));
    %y = v ./ (1 + alpha.*(1-p)./p);
end


