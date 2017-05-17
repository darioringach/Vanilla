
% Fitting the vanilla model - Spikefinder competition

global ds calcium_train spike_train X Y

for(ds=1:10)                % for each dataset
     
    dataset = num2str(ds);
    
    % read the calcium data
    calcium_train = xlsread([dataset '.train.calcium.csv']);
    calcium_train(1,:) = [];        % drop first row of cell #s
    ncells = size(calcium_train,2);
    X = nanzscore(calcium_train);   % zscore the raw data
    
    % read the spikes
    spike_train   = xlsread([dataset '.train.spikes.csv']);
    spike_train(1,:) = [];  % drop first row
    Y = spike_train;
       
    % fit the dataset -- a different case for ds==5
  

    if(ds==5)
        [x,fval] = fminsearch(@predval,[6   -.3   .8   22  -10], optimset('Display','iter'));
    else
        [x,fval] = fminsearch(@predval,[12   0   1   30 ], optimset('Display','iter'));    
    end

    
    % Print the correlation for this dataset
    
    fprintf('%s %f\n', [dataset '.train.calcium.csv'], fval);
    
    % Generate predictions
    
    spks = pred(x);
    spks = [0:ncells-1 ; spks];   % add cell numbers as first row
    dlmwrite(['pred/' dataset '.train.spikes.csv'], spks); % write the predictions
    
    % Save the optimal parameters
    save(['pred/' dataset '.train.mat'],'x','fval');  % save optimal params
    
    % Apply to test dataset if it exists
    
    if(exist([dataset '.test.calcium.csv'],'file'))   % apply to test data
        calcium_train = xlsread([dataset '.test.calcium.csv']);
        calcium_train(1,:) = []; % drop first row of cell #s
        ncells = size(calcium_train,2);
        X = nanzscore(calcium_train);
        spks = pred(x);
        spks = [0:ncells-1 ; spks];   % add cell numbers
        dlmwrite(['pred/' dataset '.test.spikes.csv'], spks);
    end
    
end
