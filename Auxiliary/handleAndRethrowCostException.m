function handleAndRethrowCostException(e, additionalCost)
% HANDLEANDRETHROWCOSTEXCEPTION  Accumulate cost in exception chain and rethrow.
%
%   HANDLEANDRETHROWCOSTEXCEPTION(E, ADDITIONALCOST) implements a cost-
%   accumulation pattern used during ODE integration failures:
%
%   1. Extracts the numeric cost stored in E.cause{end}.message (as a
%      string-encoded number, e.g. '42.5').
%   2. Adds ADDITIONALCOST to it.
%   3. Re-encodes the total as a new MException with the same identifier,
%      appends it to the cause chain, and rethrows.
%
%   This allows multiple layers of error handling to accumulate partial
%   cost contributions before bubbling up to evaluateProblem.
%   If no numeric value can be parsed from the cause message, defaults to 0.
%
%   Inputs:
%     E              - MException caught by the caller
%     ADDITIONALCOST - numeric penalty to add to the accumulated cost
%
%   See also: evaluateProblem, evaluateModel
        if ~isempty(e.cause)
            % Extract the cause message if present
            causeMessage = e.cause{end}.message;
            
            % Parse the error value from the cause message
            errorValue = str2double(causeMessage);
            if isnan(errorValue)
                errorValue = 0; % Default value if parsing fails
            end
        else
            % No cause found; set a default error value
            errorValue = 0;
        end
        
        % Add the additional value
        updatedErrorValue = errorValue + additionalCost;
        
        % Create a new message with the updated value
        newMessage = sprintf('%d', updatedErrorValue);
        
        % Create a new MException with the updated message
        newException = MException(e.identifier, newMessage);
        
        % Add the original exception as a cause
        e = e.addCause(newException);
        
        % Rethrow the exception
        throw(e);
end