function [passingFields, failingFields, errorMessages] = testRequiredFields(baseParams)
%TESTREQUIREDFIELDS Checks which fields are required for RunBakers.m
%
% Usage:
%   [pass, fail, errs] = testRequiredFields(S)
%
% Inputs:
%   S - scalar struct with all fields
%
% Outputs:
%   passingFields - cell array of field names that can be removed safely
%   failingFields - cell array of field names whose removal caused error
%   errorMessages - cell array of error messages (same length as failingFields)

    if ~isstruct(baseParams) || numel(baseParams) ~= 1
        error('Input must be a scalar struct.');
    end

    fields = fieldnames(baseParams);
    passingFields = {};
    failingFields = {};
    errorMessages = {};
    LoadData;
    baseParams.DryRun = true;
    
    for i = 1:numel(fields)
        fprintf('Testing without field "%s"...', fields{i});

        % Remove one field
        params0 = rmfield(baseParams, fields{i});
        params0.RemoveField = fields{i};

        try

            % Run the script
            clf;
            RunBakersExp;

            passingFields{end+1} = fields{i}; %#ok<AGROW>
            fprintf('Ok\n');
        catch ME
            failingFields{end+1}  = fields{i}; %#ok<AGROW>
            errorMessages{end+1}  = ME.message; %#ok<AGROW>
            fprintf('Failed with "%s"\n', ME.message);
        end
    end

    % Summary
    fprintf('\n=== Fields Causing Errors ===\n');
    for i = 1:numel(failingFields)
        fprintf('  %-20s : %s\n', failingFields{i}, errorMessages{i});
    end

    fprintf('\n=== Fields Safe to Remove ===\n');
    for i = 1:numel(passingFields)
        fprintf('  %s\n', passingFields{i});
    end
end
