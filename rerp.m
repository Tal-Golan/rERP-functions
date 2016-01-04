function [ glm_estimates, unique_event_types ] = rerp( input_data, event_latencies, event_types, window_length )
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
%

data_len = length(input_data);
unique_event_types = unique(event_types);
num_event_types = length(unique_event_types);

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

% add checks for minimal and maximal timepoints of the pulses

% build a sparse design matrix
X = sparse(i,j,1,data_len,window_length*num_event_types);

% add a constant predictor
X = cat(2,X,sparse(ones(data_len,1)));

%TODO - add artifact rejection here

% ordinary least squares: solve observedTimecourse=X*b;
b=X\input_data;

% get predicted values
glm_estimates = zeros(window_length,num_event_types);
for ii = 1:num_event_types
    glm_estimates(:,ii)=b((ii-1)*window_length+1:ii*window_length);
end

end

