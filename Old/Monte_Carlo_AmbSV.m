%% Monte Carlo Simulation of Ambiguous Decision Making
% Produces estimate of response distribution during stochastic choices in ambiguous decision making task

% Indirect Method
% Initialize variables

tic
h = waitbar(0,'Simulating...');

trialnumber = 10000;
v = 20;
p = .5;
A = .25;

Choices = zeros(trialnumber,4);

noise = .4;

for i = 1:trialnumber;                          % Monte Carlo
    alpha =  .473 + .212.*randn(1);             % Normal distribution of measured Risk Preferences
    beta = .359 + .426.*randn(1);               % Normal distribution of measured Ambiguity Preferences
    c = [0 0 0 0];
    for j = 1:4;                                % Number of trials for each value/ambiguity level condition
        if (5^alpha)*(1+noise*randn(1)) > (p*(1-beta*(A/2)))*(v^alpha)*(1+noise*randn(1));      % Subjective Value under ambiguity model
            c(j) = 0;
        else
            c(j) = 1; 
        end
    Choices(i,:) = c;
    end
    waitbar(i/trialnumber,h)
end

mChoices = mean(Choices,2);

% disp('min max mean std hist')

% min(mChoices)
% max(mChoices)
% mean(mChoices)
% std(mChoices)
hist(mChoices,50)

fin = [nnz(mChoices == 0) nnz(mChoices == .25) nnz(mChoices == .5) nnz(mChoices == .75) nnz(mChoices == 1)]./trialnumber;

disp('       0      .25      .5       .75       1')
disp(fin)
disp(['mean choice = ' num2str(mean(mChoices))])
close(h)
toc



