classdef TestFeatFloor < matlab.unittest.TestCase
    % Tests for the 'fl:' noise-floor token in evalFeatureCost.
    %
    % The floor exists because a fixed objective budget gets spent on whichever
    % term is smoothest, which is usually the one already solved
    % (Analyses/RecoveryOscillation sec. 7). It must therefore:
    %   - remove the reward below the floor, exactly, and
    %   - still defend the term above it (deadband, not deletion), and
    %   - not quietly swallow a broken evaluation (NaN).

    properties
        fd = struct('x', [10 10 10 10 10], 'v', [1 1 1 1 1]);
    end

    methods (Test)

        function noFloorIsUnchanged(t)
            fs = struct('x', [11 11 11 11 11]);
            c  = evalFeatureCost(t.fd, fs, {'x'}, 2);
            t.verifyEqual(c, 5*0.01, 'AbsTol', 1e-12);   % 5 x (1/10)^2
        end

        function floorSubtractsExactly(t)
            fs = struct('x', [11 11 11 11 11]);
            c  = evalFeatureCost(t.fd, fs, {'x|fl:0.02'}, 2);
            t.verifyEqual(c, 0.05 - 0.02, 'AbsTol', 1e-12);
        end

        function floorClampsAtZeroNotNegative(t)
            fs = struct('x', [10.1 10.1 10.1 10.1 10.1]);   % raw cost 5e-4
            c  = evalFeatureCost(t.fd, fs, {'x|fl:1.0'}, 2);
            t.verifyEqual(c, 0, 'AbsTol', 1e-12);
        end

        function aPerfectFitIsStillZero(t)
            c = evalFeatureCost(t.fd, struct('x',[10 10 10 10 10]), {'x|fl:0.02'}, 2);
            t.verifyEqual(c, 0, 'AbsTol', 1e-12);
        end

        function itIsADeadbandNotADeletion(t)
            % The whole point: below the floor costs nothing, but pushing the
            % term further out must still cost MORE. A deleted term would be
            % flat in both directions.
            near = evalFeatureCost(t.fd, struct('x',10.1*ones(1,5)), {'x|fl:0.02'}, 2);
            far  = evalFeatureCost(t.fd, struct('x',  13*ones(1,5)), {'x|fl:0.02'}, 2);
            t.verifyEqual(near, 0, 'AbsTol', 1e-12);
            t.verifyGreaterThan(far, 0.4);
        end

        function floorCombinesWithScalarWeight(t)
            fs = struct('x', [11 11 11 11 11]);
            c  = evalFeatureCost(t.fd, fs, {'x|fl:0.02|3'}, 2);
            t.verifyEqual(c, 3*(0.05 - 0.02), 'AbsTol', 1e-12);
        end

        function floorCombinesWithPerElementWeights(t)
            fs = struct('x', [11 11 11 11 11]);
            a  = evalFeatureCost(t.fd, fs, {'x|w:v'},          2);
            b  = evalFeatureCost(t.fd, fs, {'x|w:v|fl:0.02'},  2);
            t.verifyEqual(b, a - 0.02, 'AbsTol', 1e-12);
        end

        function nanPenaltySurvivesTheFloor(t)
            % A NaN is a broken evaluation, not a small residual -- it must not
            % be absorbed by the floor.
            fs = struct('x', [10 10 10 10 NaN]);
            c  = evalFeatureCost(t.fd, fs, {'x|fl:1.0'}, 2);
            t.verifyEqual(c, 1, 'AbsTol', 1e-12);          % NAN_COST = 1
        end

        function floorAppliesAfterTheSelector(t)
            % Selector first, then floor: scoring 4 segments of 5 must floor
            % against the 4 that were kept.
            fs = struct('x', [11 11 11 11 99]);
            c  = evalFeatureCost(t.fd, fs, {'x[1:4]|fl:0.02'}, 2);
            t.verifyEqual(c, 4*0.01 - 0.02, 'AbsTol', 1e-12);
        end

        function badFloorTokenErrors(t)
            t.verifyError(@() evalFeatureCost(t.fd, struct('x',ones(1,5)), ...
                {'x|fl:banana'}, 2), 'evalFeatureCost:badFloor');
        end

    end
end
