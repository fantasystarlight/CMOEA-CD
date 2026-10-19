classdef eq_DTLZ3 < PROBLEM
% <multi/many> <real> <large/none> <constrained> <expensive/none>
% Equality constrained DTLZ3
% k --- 3 --- Parameter k
% p --- 1 --- Parameter p

%------------------------------- Reference --------------------------------
% H. Jain and K. Deb, An evolutionary many-objective optimization algorithm
% using reference-point based non-dominated sorting approach, part II:
% Handling constraints and extending to an adaptive approach, IEEE
% Transactions on Evolutionary Computation, 2014, 18(4): 602-622.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2024 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    
    properties(Access = private)
        k = 3;  % Parameter k
        p = 1;  % Parameter p
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            [obj.k,obj.p] = obj.ParameterSet(3,1);
            if isempty(obj.M); obj.M = 3; end
            if isempty(obj.D); obj.D = obj.M + 9; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values and constraint violations
        function Population = Evaluation(obj,varargin)
            PopDec = varargin{1};
            PopDec = max(min(PopDec,repmat(obj.upper,size(PopDec,1),1)),repmat(obj.lower,size(PopDec,1),1));
            g      = 100*(obj.D-obj.M+1+sum((PopDec(:,obj.M:end)-0.5).^2-cos(20.*pi.*(PopDec(:,obj.M:end)-0.5)),2));
            PopObj = repmat(1+g,1,obj.M).*fliplr(cumprod([ones(size(PopDec,1),1),cos(PopDec(:,1:obj.M-1)*pi/2)],2)).*[ones(size(PopDec,1),1),sin(PopDec(:,obj.M-1:-1:1)*pi/2)];
            if obj.p == 1
                PopCon = sum((PopDec(:,1:obj.k-1)-0.5).^2, 2) - 0.4^2;
            elseif obj.p == 2
                PopCon1 = sum((PopDec(:,1:obj.k-1)-0.5).^2, 2) - 0.5^2;
                PopCon2 = sum((PopDec(:,1:obj.k-2)-0.5).^2 ,2) + (PopDec(:,obj.k-1)-1).^2 - 0.5^2;
                PopCon = [PopCon1, PopCon2];
            end
            Population = SOLUTION(PopDec,PopObj,abs(PopCon) - 1e-4,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            %R = UniformPoint(N,obj.M)/2;
            t = linspace(0,2*pi,10000);
            x = sin(t) * 0.4 + 0.5;
            y = cos(t) * 0.4 + 0.5;
            z = 0.5*ones(1,10000);
            PopDec = [x;y;z]';
            g      = 100*(3-3+1+sum((PopDec(:,3:end)-0.5).^2-cos(20.*pi.*(PopDec(:,3:end)-0.5)),2));
            R = repmat(1+g,1,3).*fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:3-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,3-1:-1:1)*pi/2)];
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if obj.M == 2
                R = obj.GetOptimum(100);
            elseif obj.M == 3
                a = linspace(0,pi/2,10)';
                R = {sin(a)*cos(a'),sin(a)*sin(a'),cos(a)*ones(size(a'))};
            else
                R = [];
            end
        end
    end
end