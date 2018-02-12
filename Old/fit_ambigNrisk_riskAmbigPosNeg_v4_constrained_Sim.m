base = 0;
    
model = 'ambigNrisk'
    
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


for s = 1:length(subjects)
    subject = subjects{s};

    filename = [subject]

    
    %ACTUAL START
    
    % get rid of no choice trials
    choice = choice_orig(find(choice_orig~=2));
    vA = vA(find(choice_orig~=2));
    AL = AL(find(choice_orig~=2));
    pA = pA(find(choice_orig~=2));


    % separate fit for pos
    choiceP = choice(find(vA>0));
    vP = vA(find(vA>0));
    pP = pA(find(vA>0));
    AP = AL(find(vA>0));
    vFP = fixed_valueP * ones(length(choiceP),1);
    pFP = fixed_prob * ones(length(choiceP),1);
    
    [info,p] = fit_ambigNrisk_model_Constrained(choiceP,vFP,vP,pFP,pP,AP,model,b0,base);

    if strcmp(model,'ambigNrisk')
        slopeP = info.b(1)
        aP = info.b(3)
        bP = info.b(2)
    elseif strcmp(model,'ambigNriskFixSlope')
        slopeP = -1;
        aP = info.b(2)
        bP = info.b(1)
    end
    r2P = info.r2
 
    % separate fir for neg
    choiceN = choice(find(vA<0));
    vN = vA(find(vA<0));
    pN = pA(find(vA<0));
    AN = AL(find(vA<0));
    vFN = fixed_valueN * ones(length(choiceN),1);
    pFN = fixed_prob * ones(length(choiceN),1);
    
    [info,p] = fit_ambigNrisk_model_Constrained(~choiceN,-vFN,-vN,pFN,pN,AN,model,b0,base);

    if strcmp(model,'ambigNrisk')
        slopeN = info.b(1)
        aN = info.b(3)
        bN = info.b(2)
    elseif strcmp(model,'ambigNriskFixSlope')
        slopeN = -1
        aN = info.b(2)
        bN = info.b(1)
    end
    r2N = info.r2

    colors =   [255 0 0;
                180 0 0;
                130 0 0;
                80 0 0;
                20 0 0; 
                52 181 233;
                7 137 247;
                3 85 155;
                ]/255;

    xP = 0:0.1:max(valueP);
    uFP = fixed_prob * (fixed_valueP).^aP;
    xN = min(valueN):0.1:0;
    uFN = -fixed_prob * (-fixed_valueN).^aN;

    all_data_subject = [valueP; riskyChoicesP; ambigChoicesP ;valueN; riskyChoicesN; ambigChoicesN];
    
    xlFile = [path 'talk_results.xls'];
    dlmwrite(xlFile, subject , '-append', 'roffset', 1, 'delimiter', ' ');  
    dlmwrite(xlFile, all_data_subject, 'coffset', 1, '-append', 'delimiter', '\t');
    
     figure
    % risk pos
    subplot(2,2,1)     
    for i = 1 :length(prob)
        plot(valueP,riskyChoicesP(i,:),'o','MarkerSize',8,'MarkerEdgeColor',colors([1 1 1])...
            ,'MarkerFaceColor',colors(i,:),'Color',colors(i,:));
          hold on
        % logistic function
        uA = prob(i) * xP.^aP;
        p = 1 ./ (1 + exp(slopeP*(uA-uFP)));

        plot(xP,p,'-','LineWidth',2,'Color',colors(i,:));
        axis([0 150 0 1])
        set(gca, 'ytick', [0 0.5 1])
    end
    title([char(subject) '  alpha gain = ' num2str(aP)]);
    
    % ambig pos
    subplot(2,2,3)
    for i = 1:length(ambig)
        plot(valueP,ambigChoicesP(i,:),'o','MarkerSize',8,'MarkerEdgeColor',colors([1 1 1]),'MarkerFaceColor',colors(length(prob)+i,:));
         hold on
% 
        % logistic function
        uA = (0.5 - bP.*ambig(i)./2) * xP.^aP;
        p = 1 ./ (1 + exp(slopeP*(uA-uFP)));


        plot(xP,p,'-','LineWidth',2,'Color',colors(length(prob)+i,:));

    end
    title([char(subject) '  beta gain = ' num2str(bP)]);

    % risk neg
    subplot(2,2,2)
        for i = 1:length(prob)
        plot(valueN,riskyChoicesN(i,:),'o','MarkerSize',8,'MarkerEdgeColor',colors([1 1 1])...
            ,'MarkerFaceColor',colors(i,:),'Color',colors(i,:));
        hold on

        % logistic function
        uA = -prob(i) * (-xN).^aN;
        p = 1 ./ (1 + exp(slopeN*(uA-uFN)));

        plot(xN,p,'-','LineWidth',2,'Color',colors(i,:));

        end
    title([char(subject) '  alpha loss = ' num2str(aN)]);
        
  
    % ambig neg
    subplot(2,2,4)
    for i = 1:length(ambig)
        plot(valueN,ambigChoicesN(i,:),'o','MarkerSize',8,'MarkerEdgeColor',colors([1 1 1]),'MarkerFaceColor',colors(length(prob)+i,:));
          hold on

        % logistic function
        uA = -(0.5 - bN.*ambig(i)./2) * (-xN).^aN;
        p = 1 ./ (1 + exp(slopeN*(uA-uFN)));


        plot(xN,p,'-','LineWidth',2,'Color',colors(length(prob)+i,:));

    end
    title([char(subject) '  beta loss = ' num2str(bN)]);

    %write into file
    fprintf(fid,'%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n',subject,aP,bP,slopeP,r2P,aN,bN,slopeN,r2N);
end

fclose(fid)
    

    

