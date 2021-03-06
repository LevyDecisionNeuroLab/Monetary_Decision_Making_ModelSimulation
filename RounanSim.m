%% Monte Carlo Simulation of Ambiguous Decision Making
% Produces estimate of response distribution during stochastic choices in ambiguous decision making task

% Indirect Method
% Initialize variables

tic
h = waitbar(0,'Running Permutation... 1');

v = [10,15,25,40,60,100];
pr = [.25,.5,.75];
a = [.2,.5,.8,1];


V = repmat(v,1,14)';
P_temp = ones(8,1)*pr;
P = [P_temp(:);.5*ones(32,1)];
A_temp = ones(8,1)*a;
A = [zeros(24,1);A_temp(:)];

%exc_index = (V==125 & P > .25 & A==0);

% V(exc_index) = [];
% P(exc_index) = [];
% A(exc_index) = [];

trialnumber = length(V);
perms = 200;
simnumbs = 100;

for reps = 1:perms

    for sims = 1:simnumbs

        Cresults = zeros(7,4);

        while sum(sum(Cresults==1)) > 12 || sum(sum(Cresults==0)) > 12

            Choices = zeros(trialnumber,1);

            alpha = .1 + 1.9.*rand(1);        %Gains Normal distribution of measured Risk Preferences
            %alpha = .97 + .75*randn(1);             %Loss
            alpha(alpha<.01) = .01;
            beta = -1 + 2.*rand(1);          %Gains Normal distribution of measured Ambiguity Preferences
            %beta = -.48 + 1.01*randn(1);            %Loss

            noise = 0+.2*randn(1);



            for i = 1:trialnumber;
                %if (10^alpha)*(1+noise*randn(1)) > (P(i)*(1-beta*(A(i)/2)))*(V(i)^alpha)*(1+noise*randn(1));      % Subjective Value under ambiguity model Gain
                if (-((10)^alpha))*(1+noise*randn(1)) > (P(i)*(1-beta*(A(i)/2)))*(-(V(i)^alpha))*(1+noise*randn(1));      % Subjective Value under ambiguity model Loss
                    Choices(i) = 0;
                else
                    Choices(i) = 1; 
                end
            end


%             Risk = Choices(1:100);
%             Ambiguity = Choices(101:160);

            %Vresults = [mean(Risk(and(V(1:100)==5,P(1:100)==.13))),mean(Risk(V(1:100)==8)),mean(Risk(V(1:100)==20)),mean(Risk(V(1:100)==50)),mean(Risk(V(1:100)==125));...
                %mean(Ambiguity(V(101:160)==5)),mean(Ambiguity(V(101:160)==8)),mean(Ambiguity(V(101:160)==20)),mean(Ambiguity(V(101:160)==50)),mean(Ambiguity(V(101:160)==125))];
            %disp(Vresults)

            Cresults = zeros(7,4);

            for valindex = 1:4
                for riskindex = 1:3
                    Cresults(riskindex,valindex) = mean(Choices(V==v(valindex) & P==pr(riskindex) & A==0));
                end
                for ambindex = 4:7
                    Cresults(ambindex,valindex) = mean(Choices(and(V==v(valindex),A==a(ambindex-3))));
                end
            end
        end


        MultiSubChoices(:,:,sims) = [V,P,A,Choices];
        MultiSubResults(:,:,sims) = Cresults;
        MultiParams(sims,:) = [alpha,beta];
        MultiNoise(sims) = noise;
        MultiFOSD(sims) = sum(Cresults(:,1))/.25;

        %waitbar(sims/simnumbs,h)
    end

    Subindex = 1:simnumbs;

    %close(h)
    %toc



    %% Fit Model

    %tic
    base = 0;

    model = 'ambigNrisk';

    fixed_valueP = 10;
    fixed_valueN = -10;

    fixed_prob = 1;

    % Initialize parameters (this could be smarter)
    if strcmp(model,'ambiguity')
       b0 = [-1 1]';
    elseif strcmp(model,'power')
       b0 = [-1 0.8];
    elseif strcmp(model, 'ambigNrisk')
       b0 = [-1 0.5 0.5]'; %slope, beta, alpha
    elseif strcmp(model, 'ambigNriskFixSlope')
       b0 = [0.5 0.5]'; %beta, alpha
    end


    for s = 1:length(Subindex)


        %ACTUAL START

        %MultiSubChoices(:,:,sims) = [V,P,A,Choices];

        % get rid of no choice trials
        choice = MultiSubChoices(:,4,s);
        %vA = MultiSubChoices(:,1,s);            %Gains
        vA = -1*MultiSubChoices(:,1,s);        %Loss
        AL = MultiSubChoices(:,3,s);
        pA = MultiSubChoices(:,2,s);

        % separate fit for pos
        choiceP = choice(find(vA>0));
        vP = vA(find(vA>0));
        pP = pA(find(vA>0));
        AP = AL(find(vA>0));
        vFP = fixed_valueP * ones(length(choiceP),1);
        pFP = fixed_prob * ones(length(choiceP),1);

        [info,p] = fit_ambigNrisk_model_Constrained(choiceP,vFP,vP,pFP,pP,AP,model,b0,base);

        if strcmp(model,'ambigNrisk')
            slopeP = info.b(1);
            aP = info.b(3);
            bP = info.b(2);
        elseif strcmp(model,'ambigNriskFixSlope')
            slopeP = -1;
            aP = info.b(2);
            bP = info.b(1);
        end
        r2P = info.r2;

        % separate fit for neg
        choiceN = choice(find(vA<0));
        vN = vA(find(vA<0));
        pN = pA(find(vA<0));
        AN = AL(find(vA<0));
        vFN = fixed_valueN * ones(length(choiceN),1);
        pFN = fixed_prob * ones(length(choiceN),1);

        [info,p] = fit_ambigNrisk_model_Constrained(~choiceN,-vFN,-vN,pFN,pN,AN,model,b0,base);

        if strcmp(model,'ambigNrisk')
            slopeN = info.b(1);
            aN = info.b(3);
            bN = info.b(2);
        elseif strcmp(model,'ambigNriskFixSlope')
            slopeN = -1;
            aN = info.b(2);
            bN = info.b(1);
        end
        r2N = info.r2;

        %ParamsResults(s,:) = [aP,bP,slopeP,r2P];
        ParamsResults(s,:) = [aN,bN,slopeN,r2N];

        %waitbar(s/length(Subindex),h,'Fitting...')

    end
    
    Choice4d(:,:,:,reps) = MultiSubChoices;
    SubResults4d(:,:,:,reps) = MultiSubResults;
    gParams4d(:,:,reps) = MultiParams;
    Noise4d(:,reps) = MultiNoise;
    fosd4d(:,reps) = MultiFOSD;
    fitParams4d(:,:,reps) = ParamsResults;
    
    waitbar(reps/perms,h,['Runing Permutation... ' num2str(reps+1)])
end


close(h)    
toc




