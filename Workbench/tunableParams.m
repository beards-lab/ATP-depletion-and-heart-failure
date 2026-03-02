function flattenedOutput = tunableParams(params0)
% Parameters to exclude
excludeParams = {
 'baseline_dummy',
 'dS',
'A1AttachmentWidth', 
'Amax', 
'BreakOnODEUnstable', 
'Ca', 
'DryRun', 
'EvalAtp', 
'EvalFeatures', 
'EvalFitSlackOnset', 
'FV_velocities', 
'F_act_UseP31', 
'FudgeVmax', 
'LSE0', 
'LXBpivot', 
'LegacyStrainFlipping', 
'Lsc0', 
'ML', 
'MatchTimeSegments', 
'MaxRunTime', 
'MaxSlackNegativeForce', 
'MaxSpaceExtensionCount', 
'MaxStrainArraySize', 
'MgADP', 
'MgATP', 
'N', 
'NumberOfStates', 
'OptimizeFVInit', 
'OutputAtSL', 
'PU0', 
'Pi', 
'PlotEachSeparately', 
'PlotFullscreen', 
'PlotProbsOnFig', 
'ReduceSpace', 
'RescaleOutputDt', 
'RunForceLengthEstim', 
'RunForceVelocity', 
'RunForceVelocityTime', 
'RunKtr', 
'RunSlack', 
'RunSlackSegments', 
'RunStairs', 
'SL0', 
'SLmax', 
'SaveBest', 
'ShowArrayShiftWarnings', 
'ShowResidualPlots', 
'ShowStatePlots', 
'SimTitle', 
'Slim', 
'Slim_l', 
'Slim_r', 
'StrainExp', 
'UseA1AttachmentKernel', 
'UseA2MechanicalRecocking', 
'UseA2Popping', 
'UseA2Reattaching', 
'UseAtpOnUNR', 
'UseCa', 
'UseCalculatedN', 
'UseDirectSRXTransition', 
'UseForceOnsetShift', 
'UseKstiff3', 
'UseKtrProtocol', 
'UseLimitedSRDTransition', 
'UseMutualPairingAttachment', 
'UseNegativeForceRip', 
'UseNegativeKstiff', 
'UseOverlap', 
'UseOverlapFactor', 
'UseP31Shift', 
'UsePassive', 
'UsePassiveForSR', 
'UsePieceWiseStrainDep', 
'UseSLInput', 
'UseSafeSRReset', 
'UseSerialStiffness', 
'UseSlack', 
'UseSpaceDiscretization', 
'UseSpaceExtension', 
'UseSpaceInterpolation', 
'UseStrainDep4R1D', 
'UseStrictDetachmentAt', 
'UseSuperRelaxed', 
'UseSuperRelaxedADP', 
'UseTORNegShift', 
'UseTitinIdentifiedPassive', 
'UseTitinInterpolation', 
'UseTitinModel', 
'UseUniformTransitionFunc', 
'UseWDetachment', 
'ValuesInTime', 
'Velocity', 
'WindowsOverflowStepCount', 
'datatable', 
'drawForceOnset', 
'g', 
's',
'fn',
'gamma', 
'ghostLoad', 
'ghostSave', 
'justPlotStateTransitionsFlag', 
'modelFcn', 
'mods', 
'ss', 
'useHalfActiveForSR', 
'velocitytableonfile', 
'PieceWiseStrainDepParams__2'   ,
'PieceWiseStrainDepParams__3'   ,
'PieceWiseStrainDepX__2'        ,
'PieceWiseStrainDepX__3'        ,
'PieceWiseStrainDep2X__2'       ,
'PieceWiseStrainDep2X__3'       ,
'PieceWiseStrainDep2X__4'       ,
'PieceWiseStrainDep2Params__2'  ,
'PieceWiseStrainDep2Params__3'  ,
'PieceWiseStrainDep2Params__4'  ,
'PieceWiseStrainDepR1DX__2'     ,
'PieceWiseStrainDepR1DX__3'     ,
'PieceWiseStrainDepR1DParams__2',
'PieceWiseStrainDepR1DParams__3',
'PieceWiseStrainDepR21X__2'     ,
'PieceWiseStrainDepR21X__3'     ,
'PieceWiseStrainDepR21Params__2',
'PieceWiseStrainDepR21Params__3',
'vmax',
'UseMaxwellDashpot'
};

% Initialize the output struct
flattenedOutput = struct();

fields = fieldnames(params0);

for i = 1:numel(fields)
    fieldName = fields{i};
    val = params0.(fieldName);
    
    % 1. Check if the field is in the exclusion list
    if ismember(fieldName, excludeParams)
        continue;
    end
    
    % 2. Check if it's a single numerical value or a boolean
    if (isnumeric(val) || islogical(val)) && isscalar(val)
        flattenedOutput.(fieldName) = val;        
    % 3. Check if it's a numerical array (excluding strings/cells)
    elseif isnumeric(val) && ~isscalar(val) && isvector(val)
        for idx = 1:length(val)
            indexedName = sprintf('%s__%d', fieldName, idx);
                if ismember(indexedName, excludeParams)
                    continue;
                end
            flattenedOutput.(indexedName) = val(idx);
        end
    else
        fprintf('Skipping %s.. \n', fieldName)
    end
    
    % All other types (structs, strings, etc.) are implicitly ignored
end

% Display Result
% disp(flattenedOutput);
end