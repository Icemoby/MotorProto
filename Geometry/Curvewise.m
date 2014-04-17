classdef Curvewise < Geometry2D
    properties (SetAccess=protected)
      	curves
        orientation
        domains
    end
    methods
        function obj=Curvewise(varargin)
            obj = obj@Geometry2D;
            n = numel(varargin);
            for i = 1:2:n
                obj.(varargin{i})=varargin{i+1};
            end
            n = numel(obj.curves);
            for i = 1:n
                obj.curves{i} = obj.curves{i}.clone;
            end
            obj.refresh;
        end
        
        function refresh(obj)
            obj.sortCurves;
            n = numel(obj.curves);
            for i=1:n
                obj.curves{i}.refresh
            end
            notify(obj,'Geometry2DUpdate');
        end
        
        function build(~)
        end
        
      	function [In, On, N] = inOn(objIn,X,Y)
            c = objIn.curves;
            n = numel(c);
            m = length(X);
            InC = false(m,n);
            OnC = false(m,n);
            NC = zeros(m,2*n);
            for i = 1:n
                [InC(:,i), OnC(:,i), NC(:,2*i-1:2*i)] = c{i}.inOn(X,Y);
            end
            In = all(InC,2);
            On = false(m,1);
            N = zeros(m,2);
            for i = 1:n
                OnI = (all(InC(:,[1:i-1,i+1:n]),2)&OnC(:,i));
                N(OnI,:) = NC(OnI,2*i-1:2*i);
                On = On|OnI;
            end
        end
        
        function objOut = clone(objIn,cType)
            if nargin==1
                cType = 'parameter';
            end
            objOut = curvewise;
            c = objIn.curves;
            n = numel(c);
            for i=1:n
                c{i} = c{i}.clone(cType);
            end
            objOut.curves = c;
        end
    end
end