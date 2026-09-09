classdef TestFeatSelector < matlab.unittest.TestCase
    % TESTFEATSELECTOR  The [selector] suffix on a feature name.
    %
    % Segment ranges ('rsK80[1:4]') exist so a slack cycle that is a different
    % manoeuvre can be dropped from a feature without dropping the feature --
    % see Analyses/RecoveryOscillation. The list form is the trap that made
    % this worth a test: str2double('1,3,5') is 135, because MATLAB reads the
    % comma as a thousands separator, so a list must be matched BEFORE any
    % numeric parse or it silently becomes an out-of-range index.

    properties (Constant)
        v = [10 20 30 40 50];
    end

    methods (Test)
        function emptySelectorPassesThrough(t)
            t.verifyEqual(applyFeatSelector(t.v, ''), t.v);
        end

        function singleIndex(t)
            t.verifyEqual(applyFeatSelector(t.v, '1'), 10);
            t.verifyEqual(applyFeatSelector(t.v, '5'), 50);
        end

        function outOfRangeSingleIsNaN(t)
            t.verifyTrue(isnan(applyFeatSelector(t.v, '9')));
        end

        function contiguousRange(t)
            t.verifyEqual(applyFeatSelector(t.v, '1:4'), [10 20 30 40]);
            t.verifyEqual(applyFeatSelector(t.v, '2:4'), [20 30 40]);
        end

        function rangeIsClippedNotAnError(t)
            t.verifyEqual(applyFeatSelector(t.v, '4:9'), [40 50]);
            t.verifyTrue(isnan(applyFeatSelector(t.v, '7:9')));
        end

        function listIsNotReadAsAThousandsSeparator(t)
            % str2double('1,3,5') == 135. If the list form is parsed after
            % the scalar form this returns NaN instead of three segments.
            t.verifyEqual(applyFeatSelector(t.v, '1,3,5'), [10 30 50]);
            t.verifyEqual(applyFeatSelector(t.v, '1:3,5'), [10 20 30 50]);
        end

        function reductionsStillWork(t)
            t.verifyEqual(applyFeatSelector(t.v, 'end'),    50);
            t.verifyEqual(applyFeatSelector(t.v, 'mean'),   30);
            t.verifyEqual(applyFeatSelector(t.v, 'median'), 30);
            t.verifyEqual(applyFeatSelector(t.v, 'min'),    10);
            t.verifyEqual(applyFeatSelector(t.v, 'max'),    50);
        end

        function parseSplitsNameFromSelector(t)
            [n, s] = parseFeatSelector('rsK80[1:4]');
            t.verifyEqual(n, 'rsK80'); t.verifyEqual(s, '1:4');
            [n2, s2] = parseFeatSelector('rsK80');
            t.verifyEqual(n2, 'rsK80'); t.verifyEqual(s2, '');
        end

        function bracketCommaIsNotAGroupSeparator(t)
            % A grouped fn entry is split on top-level commas. A comma inside
            % a segment list must not split the entry into two terms.
            t.verifyFalse(hasTopLevelComma('rsK80[1,3,5]|0.6'));
            t.verifyTrue( hasTopLevelComma('k2|100-1500|0.01,ksrd|0-20|0.001'));
            p = splitTopLevelComma('rsK80[1,3,5]|0.6');
            t.verifyEqual(numel(p), 1);
        end

        function costUsesOnlyTheSelectedSegments(t)
            % The whole point: a wrong segment 5 must not reach the cost.
            fd.rsK80 = [50 50 50 50 50];
            fs.rsK80 = [50 50 50 50 500];
            cAll = evalFeatureCost(fd, fs, {'rsK80'},      2);
            c14  = evalFeatureCost(fd, fs, {'rsK80[1:4]'}, 2);
            t.verifyGreaterThan(cAll, 1);
            t.verifyEqual(c14, 0, 'AbsTol', 1e-12);
        end
    end
end
