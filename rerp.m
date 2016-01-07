function [ glm_estimates, design_matrix ] = rerp( input_data, varargin)
%RERP performs a linear regression based analysis of event-related data.
% 
% Usage: 
% b = RERP(data, latencies, types, len) where data is the observed time 
% course, latencies is a vector holding the event onset latencies, types is
% a vector of the same length where each value codes the type ("condition")
% of the corresponding event, and len is the length of the window to be 
% analyzed, returns the matrix b where each column is the beta coefficients 
% associated with the timecourse of a specific event type. 
% len can either be a scalar for a fixed analysis window length, or a
% vector with the same number of elements as there are unique event types
% to indicate a different window length for each event type (order is
% expected to correspond to ascending event type values). In the case of
% variable window lengths, columns of the b matrix will be padded with NaNs
% for analysis windows shorter than the maximum length. 
% b = RERP(data, latencies, types, len, art) where art is a vector of time
% points which should be excluded from analysis (artifacts), will compute 
% the regression after removing the corresponding data points from the data
% and the design matrix. 
% [b, D] = RERP(...) will also return the struct D holding the resulting
% design matrix and meta-data used in the regression analysis. 
% b = RERP(data, D) where is the design matrix returned by a previous run
% of the function, will apply the given design matrix to the new data. 
% b = RERP(data, D, art) will also remove from analysis the data points 
% indicated by the vector art. 
% 
% For theoretical background, see: 
% 1. Smith, N.J. & Kutas, M., 2014. Regression-based estimation of ERP 
%    waveforms: I. The rERP framework. Psychophysiology, 52, pp.157–168.
% 2. Smith, N.J. & Kutas, M., 2015. Regression-based estimation of ERP 
%    waveforms: II. Nonlinear effects, overlap correction, and practical 
%    considerations. Psychophysiology, 52, pp.169–181.
% 
% Written by Tal Golan, Edden Gerber and Tamar Regev, Jan. 2016 

tictoc = tic;

% Handle input 
nargin = length(varargin);

% If additional predictors exist assign them to variable
if nargin>1 && length(input_data)==length(varargin{nargin})
    extra_predictors=varargin{nargin};
    nargin=nargin-1; % ignore last input in counting input variables
else 
    extra_predictors=[];
end

if nargin < 1
    error('Not enough input arguments');
    
elseif nargin < 3 % 1 or 2 additional arguments
    x_exists = true;
    design_matrix = varargin{1};
    event_latencies = design_matrix.event_latencies;
    event_types = design_matrix.event_types;
    window_length = design_matrix.window_length;
    artifact_indexes = [];
    if nargin == 2 % artifact indexes included
        artifact_indexes = varargin{2};
    end
    
elseif nargin < 5 % 3 or 4 additional arguments
    x_exists = false;
    event_latencies = varargin{1};
    event_types = varargin{2};
    window_length = varargin{3};
    artifact_indexes = []; 
    if nargin == 4 % artifact indexes included
        artifact_indexes = varargin{4};
    end
    
else
    error('Too many input arguments');
end

% Initialize variables
data_len = length(input_data);
unique_event_types = unique(event_types);
num_event_types = length(unique_event_types);
% In case window_length is a scalar, duplicate it for all event types
if isscalar(window_length)
    window_length = repmat(window_length, num_event_types, 1);
end
% In case unique_event_types are not consecutive integers
event_type_indexes = zeros(size(event_types));
for ii = 1:num_event_types 
    event_type_indexes(event_types == unique_event_types(ii)) = ii;
end
% Make sure all vectors are columns
if isrow(event_latencies); event_latencies = event_latencies'; end
if isrow(event_types); event_types = event_types'; end
if isrow(window_length); window_length = window_length'; end

% Check if any event windows exceed data length. Warn and fix by eliminating them :
events_exceeding = event_latencies + window_length(event_type_indexes) >= data_len;
if any(events_exceeding)
    warning([num2str(sum(events_exceeding)) ' events exceed data length'])
    warning('eliminating events from analysis')
    event_latencies(events_exceeding)=[];
    event_types(events_exceeding)=[];
end

% Calculate design matrix, if not provided as input
if x_exists
    disp('Design matrix provided as input');
    X = design_matrix.matrix;
else
    disp('Calculating design matrix...');
    
    i = nan(sum(window_length(event_type_indexes)),1);
    j = nan(sum(window_length(event_type_indexes)),1);

    n = 1;
    for ii = 1:num_event_types
        current_event_type = unique_event_types(ii);
        num_events = sum(event_types==current_event_type);
        
        pulse_latencies=bsxfun(@plus,event_latencies(event_types==current_event_type),0:window_length(ii)-1); % (nEvents x nLags)
        pulse_predictors=repmat(1:window_length(ii),numel(event_latencies(event_types==current_event_type)),1); % (nEvents x nLags)
        
        i(n:n+num_events*window_length(ii)-1) = pulse_latencies(:);
        j(n:n+num_events*window_length(ii)-1) = pulse_predictors(:) + sum(window_length(1:ii-1));
        
        n = n + num_events*window_length(ii);
    end
    % Build a sparse design matrix
    X = sparse(i,j,1,data_len,sum(window_length));
    % Add a constant predictor
    X = cat(2,X,sparse(ones(data_len,1)));
    % Add extra predictors
    if length(extra_predictors)==size(extra_predictors,1)
        X=cat(2,X,sparse(extra_predictors));
        num_extra_predictors=size(extra_predictors,2);
    elseif length(extra_predictors)==size(extra_predictors,2)
        X=cat(2,X,sparse(extra_predictors'));
        num_extra_predictors=size(extra_predictors,1);
    else
        warning('Wrong dimensions of additional predictors matrix');
    end
    
    % Assign design matrix struct
    design_matrix = struct;
    design_matrix.matrix = X;
    design_matrix.event_latencies = event_latencies;
    design_matrix.event_types = event_types;
    design_matrix.unique_event_types = unique_event_types;
    design_matrix.num_event_types = num_event_types;
    design_matrix.window_length = window_length;
    design_matrix.num_extra_predictors=num_extra_predictors;
end

% Remove artifacts from analysis
X(artifact_indexes,:) = [];
input_data(artifact_indexes,:) = [];

% Ordinary least squares: solve observedTimecourse=X*b;
disp('running regression...');
b = X\input_data;

% Parse beta vector into segments
glm_estimates = nan(max(window_length),num_event_types+min(size(extra_predictors)));
for ii = 1:num_event_types
    glm_estimates(1:window_length(ii),ii) = b(sum(window_length(1:ii-1))+1:sum(window_length(1:ii)));
%     b((ii-1)*window_length+1:ii*window_length);
end
for ii=1:min(size(extra_predictors)) 
   glm_estimates(1,num_event_types+ii) = b(num_event_types*(window_length(1))+ii);
end

disp(['Done in ' num2str(toc(tictoc)) ' sec.']);

end

