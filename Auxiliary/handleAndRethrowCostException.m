function handleAndRethrowCostException(e, additionalCost)
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