% CHOICE_PROB_AMBIGNRISK                Binary logit choice probability
% 
%     p = choice_prob_ambigNrisk(vF,vA,pF,pA,AL,beta,model);
%
%     INPUTS
%     vF    - values of fixed option
%     vA    - values of ambiguous option
%     pF    - probability of fixed option
%     pA    - probability of ambiguous option
%     AL    - ambiguity level
%     beta  - Parameters corresponding to MODEL
%     model - String indicating which model to fit; currently valid are:
%               'ambigNrisk' - (p-beta(2)*AL/2)*v^beta(3*)
%
%     OUTPUTS
%     p     - choice probabilities for the *SHORTER* option
%

%
%     REVISION HISTORY:
%     brian 03.14.06 written
%     ifat  12.01.06 adapted for ambiguity and risk

function p = choice_prob_ambigNrisk(base,vF,vA,pF,pA,AL,beta,model);

if strcmp(model,'ambigNriskPosNeg2slopes')
    % use 2 logistics for pos and neg
    % change that!!!
%     vFp = zeros(size(vF));
%     vFn = zeros(size(vF));
%     vAp = zeros(size(vA));
%     vAn = zeros(size(vA));
%     vFp(find(vF>0)) = vF(find(vF>0));
%     vFn(find(vF<0)) = vF(find(vF<0));
%     vAp(find(vA>0)) = vA(find(vA>0));
%     vAn(find(vA<0)) = vA(find(vA<0));
% 
%     uF = ambig_utility_POSNEG(vFp,vFn,pF,zeros(size(pF)),beta(3),beta(2),beta(4),beta(5),beta(6),model);
%     uA = ambig_utility_POSNEG(vAp,vAn,pA,AL,beta(3),beta(2),beta(4),beta(5),beta(6),model);
%     slopeP = beta(1);
%     slopeN = beta(7);
%     
%     p1 = 1 ./ (1 + exp(beta(1)*(uA-uF)));
else
    % single logistic
    if (strcmp(model,'ambigNriskFixSlope'))
        uF = ambig_utility(vF,pF,zeros(size(vF)),beta(2),beta(1),model); %fixed non-ambiguous
        uA = ambig_utility(vA,pA,AL,beta(2),beta(1),model); % ambiguous
        slope = -1
    elseif (strcmp(model,'ambigNrisk') || strcmp(model,'ambiguity') || strcmp(model,'ambigPower'))
        % this is the one we are currently using!
        uF = ambig_utility(base,vF,pF,zeros(size(vF)),beta(3),beta(2),model); %fixed non-ambiguous
        uA = ambig_utility(base,vA,pA,AL,beta(3),beta(2),model); % ambiguous
        slope = beta(1);
    elseif strcmp(model,'power')
        uF = utility(vF,pF,beta(2),model); % fixed
        uA = utility(vA,pA,beta(2),model); % risky
        slope = beta(1);
    elseif strcmp(model,'ambigWeighted')
        uF = ambigWeighted_utility(vF,pF,zeros(size(vF)),beta(3),beta(2),beta(4),model);
        uA = ambigWeighted_utility(vA,pA,AL,beta(3),beta(2),beta(4),model);
        slope = beta(1);
    elseif  strcmp(model,'discounting')
        uF = discounting_utility(vF,pF,zeros(size(vF)),beta(3),beta(2),beta(4),model); %fixed non-ambiguous
        uA = discounting_utility(vA,pA,AL,beta(3),beta(2),beta(4),model); % ambiguous
        slope = beta(1);
    elseif (strcmp(model,'ambigNriskPosNeg')) 
        % fit will be done simultaneously for pos and neg but with one slope
        vFp = zeros(size(vF));
        vFn = zeros(size(vF));
        vAp = zeros(size(vA));
        vAn = zeros(size(vA));
        vFp(find(vF+base>0)) = vF(find(vF+base>0));
        vFn(find(vF+base<0)) = vF(find(vF+base<0));
        vAp(find(vA+base>0)) = vA(find(vA+base>0));
        vAn(find(vA+base<0)) = vA(find(vA+base<0));
        uF = ambig_utility_POSNEG(base,vFp,vFn,pF,zeros(size(pF)),beta(3),beta(2),beta(5),beta(4),beta(6),model);
        uA = ambig_utility_POSNEG(base,vAp,vAn,pA,AL,beta(3),beta(2),beta(5),beta(4),beta(6),model);
        slope = beta(1);
    elseif (strcmp(model,'ambigNriskPosNegOneAlpha')) 
        % fit will be done simultaneously for pos and neg but with one slope
        vFp = zeros(size(vF));
        vFn = zeros(size(vF));
        vAp = zeros(size(vA));
        vAn = zeros(size(vA));
        vFp(find(vF>0)) = vF(find(vF>0));
        vFn(find(vF<0)) = vF(find(vF<0));
        vAp(find(vA>0)) = vA(find(vA>0));
        vAn(find(vA<0)) = vA(find(vA<0));
        uF = ambig_utility_POSNEG(base,vFp,vFn,pF,zeros(size(pF)),beta(3),beta(2),beta(3),beta(4),beta(5),model);
        uA = ambig_utility_POSNEG(base,vAp,vAn,pA,AL,beta(3),beta(2),beta(3),beta(4),beta(5),model);
        slope = beta(1);
    
    end
    %s = ones(size(uA)); %sign(uA);
    p = 1 ./ (1 + exp(slope*(uA-uF)));
end

return

