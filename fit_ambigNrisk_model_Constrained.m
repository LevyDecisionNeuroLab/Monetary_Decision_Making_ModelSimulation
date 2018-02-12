% FIT_AMBIGNRISK_MODEL      Fit a variety of probabilistic ambiguity models
% 
%     [info,p] = fit_ambigNrisk_model(choice,vF,vA,pF,pA,AL,model);
%
%     Fits a binary logit model by maximum likelihood.
%
%     INPUTS
%     choice      - Dependent variable. The data should be *ungrouped*,
%                   such that CHOICE is a column of 0s and 1s, where 1 indicates 
%                   a choice of the FIXED option.
%     vF          - value of fixed lottery
%     vA          - value of ambiguous lottery
%     pF          - probability of fixed lottery
%     pA          - probability of ambiguous lottery
%     AL          - ambiguity level
%     model       - String indicating which model to fit; currently valid are:
%                   'ambigNrisk'   - power with subjective probability, estimates both risk and ambiguity coefficients 
%                    
%                   Multiple models can be fit by passing in a cell array
%                   of strings. 
%
%     OUTPUTS
%     info       - data structure with following fields:
%                     .nobs      - number of observations
%                     .nb        - number of parameters
%                     .optimizer - function minimizer used
%                     .exitflag  - see FMINSEARCH
%                     .b         - fitted parameters; note that for all the
%                                  available models, the first element of B
%                                  is a noise term for the logistic
%                                  function, the remaining elements are
%                                  parameters for the selected discount
%                                  functions. eg., for model='exp', B(2) is
%                                  the time constant of the exponential
%                                  decay.
%                     .LL        - log-likelihood evaluated at maximum
%                     .LL0       - restricted (minimal model) log-likelihood
%                     .AIC       - Akaike's Information Criterion 
%                     .BIC       - Schwartz's Bayesian Information Criterion 
%                     .r2        - pseudo r-squared
%                   This is a struct array if multiple models are fit.
%     p           - Estimated choice probabilities evaluated at the values
%                   delays specified by the inputs vS, vR, dS, dL. This is
%                   a cell array if multiple models are fit.
%
%     EXAMPLES
%     see TEST_FAKE_DATA_AMBIGUITTY, TEST_FAKE_DATA, TEST_JOE_DATA, and TEST_KENWAY_DATA

%
%     REVISION HISTORY:
%     brian 03.10.06 written
%     brian 03.14.06 added fallback to FMINSEARCH, multiple fit capability
%     ifat  12.01.06 adapted for ambiguity and risk + CI

function [info,p] = fit_ambigNrisk_model(choice,vF,vA,pF,pA,AL,model,b0,base);
% If multiple model fits requested, loop and pack

if iscell(model)
   for i = 1:length(model)
      [info(i),p{i}] = fit_ambigNrisk_model(choice,vF,vA,pF,pA,AL,model{i});
   end
   return;
end

thresh = 0.05;
nobs = length(choice);

% Fit model, attempting to use FMINUNC first, then falling back to FMINSEARCH
% if exist('fminunc','file')
%    try
%       optimizer = 'fminunc';
%       OPTIONS = optimset('Display','off','LargeScale','off');
%           [b,negLL,exitflag,convg,g,H] = fminunc(@local_negLL,b0,OPTIONS,choice,vF,vA,pF,pA,AL,model,base);
%       
%       if exitflag ~= 1 % trap occasional linesearch failures
%          optimizer = 'fminsearch';
%          fprintf('FMINUNC failed to converge, switching to FMINSEARCH\n');
%       end         
%    catch
%       optimizer = 'fminsearch';
%       fprintf('Problem using FMINUNC, switching to FMINSEARCH\n');
%    end
% else
%    optimizer = 'fminsearch';
% end
% 
% if strcmp(optimizer,'fminsearch')
%    optimizer = 'fminsearch';
%    OPTIONS = optimset('Display','off','TolCon',1e-6,'TolFun',1e-5,'TolX',1e-5,...
%       'DiffMinChange',1e-4,'Maxiter',100000,'MaxFunEvals',20000);
%    [b,negLL,exitflag,convg] = fminsearch(@local_negLL,b0,OPTIONS,choice,vF,vA,pF,pA,AL,model,base);
% end
% 
% if exitflag ~= 1
%    fprintf('Optimization FAILED, #iterations = %g\n',convg.iterations);
% else
%    fprintf('Optimization CONVERGED, #iterations = %g\n',convg.iterations);
optimizer = 'fmincon';

OPTIONS = optimset('Display','off','TolCon',1e-6,'TolFun',1e-5,'TolX',1e-5,...
    'DiffMinChange',1e-4,'Maxiter',100000,'MaxFunEvals',20000,'Algorithm','active-set');
[b,negLL,exitflag,convg] = fmincon(@local_negLL,b0,[],[],[],[],[-inf -3.67 .0894],[inf 4 4.34],[],OPTIONS,choice,vF,vA,pF,pA,AL,model,base);
%X = fmincon(fun,x0,A,B,Aeq,Beq,lb,ub,nonlcon,Opts,a,b)



% Choice probabilities (for VARIED)
p = choice_prob_ambigNrisk(base,vF,vA,pF,pA,AL,b,model);

% Unrestricted log-likelihood
LL = -negLL;
% Restricted log-likelihood
LL0 = sum((choice==1).*log(0.5) + (1 - (choice==1)).*log(0.5));

% Confidence interval, requires Hessian from FMINUNC
try
    invH = inv(-H);
    se = sqrt(diag(-invH));
catch
end

info.nobs = nobs;
info.nb = length(b);
info.model = model;
info.optimizer = optimizer;
info.exitflag = exitflag;
info.b = b;

try
    info.se = se;
    info.ci = [b-se*norminv(1-thresh/2) b+se*norminv(1-thresh/2)]; % Wald confidence
    info.tstat = b./se;
catch
end

info.LL = LL;
info.LL0 = LL0;
info.AIC = -2*LL + 2*length(b);
info.BIC = -2*LL + length(b)*log(nobs);
info.r2 = 1 - LL/LL0;

%----- LOCAL FUNCTIONS
function sumerr = local_negLL(beta,choice,vF,vA,pF,pA,AL,model,base);
p = choice_prob_ambigNrisk(base,vF,vA,pF,pA,AL,beta,model);

% Trap log(0)
ind = p == 1;
p(ind) = 0.9999;
ind = p == 0;
p(ind) = 0.0001;
% Log-likelihood
err = (choice==1).*log(p) + (1 - (choice==1)).*log(1-p);
% Sum of -log-likelihood
sumerr = -sum(err);



