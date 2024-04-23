classdef frac_parameters
    properties
        a
        d
        Tf
        N
    end

    methods
        function obj = frac_parameters(alpha, delta, Tf, N)
            obj.a = alpha;
            obj.d = delta;
            obj.Tf = Tf;
            obj.N = N;
        end
    end
end