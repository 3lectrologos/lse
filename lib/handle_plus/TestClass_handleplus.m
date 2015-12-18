classdef TestClass_handleplus < handleplus
    %TESTCLASS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        prop1 = 1;
        prop2 = 2;
    end
    properties (Dependent)
        prop3
    end
    
    methods
        function IncrementVersion(obj)
            obj.version = obj.version+1;
        end
    end
    
    
end

