function [ glm_estimates, design_matrix ] = rerp( input_data, varargin)
%
%RERP performs a linear regression based analysis of event-related data.
% 
% Usage: 
% b = RERP(data, latencies, types, len) where data is the observed time 
% course, latencies is a vector holding the event onset latencies, types is
% a vector of the same length where each value codes the type ("condition")
% of the corresponding event, and len is the length of the window to be 
% analyzed, returns the matrix b where each column is the beta coefficients 
% associated with the timecourse of a specific event type. 
% [b, t] = RERP(..) also returns t, a list of the unique event types, where
% event type t(i) corresponds to the timecourse b(:,i). 
% 
% For theoretical background, see: 
% 1. Smith, N.J. & Kutas, M., 2014. Regression-based estimation of ERP 
%    waveforms: I. The rERP framework. Psychophysiology, 52, pp.157–168.
% 2. Smith, N.J. & Kutas, M., 2015. Regression-based estimation of ERP 
%    waveforms: II. Nonlinear effects, overlap correction, and practical 
%    considerations. Psychophysiology, 52, pp.169–181.
% 
% Written by Tal Golan and Edden Gerber, Jan. 2016 
% - Tamar Regev - 5/1/2016 
%   - added warning and fix if last event window exceeds data length
%   - added artifact rejection 
%   - optionally accept matrix X as input and do not calculate it 
%     --> input ossibilities are
%       rerp( input_data, design_matrix)
%       rerp( input_data, design_matrix, artifact_indexes)
%       rerp( input_data, event_latencies, event_types, window_length)
%       rerp( input_data, event_latencies, event_types, window_length, artifact_indexes )

tictot = tic;

switch length(varargin)
    case 0
        error('not enough input arguments')
    case 1
        Xexists = true;
        design_matrix = varargin{1};
        event_types = design_matrix.event_types;
        event_latencies = design_matrix.event_latencies;
        window_length = design_matrix.window_length;
    case 2
        Xexists = true;
        design_matrix = varargin{1};
        artifact_indexes = varargin{2};
        event_types = design_matrix.event_types; 
        event_latencies = design_matrix.event_latencies;
        window_length = design_matrix.window_length;
    case 3
        Xexists = false;
        event_latencies = varargin{1};
        event_types = varargin{2};
        window_length = varargin{3};
    case 4
        Xexists = false;
        event_latencies = varargin{1};
        event_types = varargin{2};
        window_length = varargin{3};
        artifact_indexes = varargin{4};
end

data_len = length(input_data);
unique_event_types = unique(event_types);
num_event_types = length(unique_event_types);

% Check if last event windows exceed data length. Warn and fix by eliminating it :
nEvents_exceeding = sum(event_latencies+window_length-1 >= data_len);
if nEvents_exceeding
    warning([num2str(nEvents_exceeding) 'last events exceed data length'])
    warning('eliminating events from analysis')
    event_latencies(end+1-nEvents_exceeding:end)=[];
end

%calculate design matrix, if not provided as input
if Xexists
    disp('Design matrix provided as input')
    X = design_matrix.matrix;
    window_length = design_matrix.window_length;
    num_event_types = design_matrix.num_event_types;
else
    disp(['Calculating design matrix...'])
    
    i = nan(length(event_latencies)*window_length,1);
    j = nan(length(event_latencies)*window_length,1);

    n = 1;
    for ii = 1:num_event_types
        current_event_type = unique_event_types(ii);
        num_events = sum(event_types==current_event_type);

        pulse_latencies=bsxfun(@plus,event_latencies(event_types==current_event_type),0:window_length-1); % (nEvents x nLags)
        pulse_predictors=repmat(1:window_length,numel(event_latencies(event_types==current_event_type)),1); % (nEvents x nLags)

        i(n:n+num_events*window_length-1) = pulse_latencies(:);
        j(n:n+num_events*window_length-1) = pulse_predictors(:) + window_length*(ii-1);

        n = n + num_events*window_length;

    end
    
    % build a sparse design matrix
    X = sparse(i,j,1,data_len,window_length*num_event_types);

    % add a constant predictor
    X = cat(2,X,sparse(ones(data_len,1)));
    disp(['Done'])
    
    %assign design matrix struct
    design_matrix = struct;
    design_matrix.matrix = X;
    design_matrix.event_types = event_types;
    design_matrix.unique_event_types = unique_event_types;
    design_matrix.window_length = window_length;
    design_matrix.num_event_types = num_event_types;
    design_matrix.event_latencies = event_latencies;
end

%reject artifacts
X(artifact_indexes,:) = [];
input_data(artifact_indexes,:) = [];

% ordinary least squares: solve observedTimecourse=X*b;
disp('running regression...');tic;
b = X\input_data;

% get predicted values
glm_estimates = zeros(window_length,num_event_types);
for ii = 1:num_event_types
    glm_estimates(:,ii)=b((ii-1)*window_length+1:ii*window_length);
end
disp(['Done in ' num2str(toc) ' sec.'])
disp(['Done all in ' num2str(toc(tictot)) ' sec.'])
end

