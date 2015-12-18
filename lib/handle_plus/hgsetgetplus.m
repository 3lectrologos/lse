classdef hgsetgetplus < handleplus & hgsetget
% 
% This class is an extension to the hgsetget class, which allows for deep copies of objects. If OBJ is a member of
% "hgsetgetplus", a deep copy may be created using:
%
% >> OBJ_new = OBJ.copy;
% >> OBJ_new = OBJ.clone;  % Alias-method, same as OBJ.copy;
%
% In the call to the copy-method, valid properties may be specified to override existing values of the copied object:
% 
% OBJ = My_hgsetgetplus_class('prop1', 1,'prop2','hallo'); % Assuming prop1 & prop2 to be valid properties of the class
% OBJ_new = OBJ.copy('prop1',123);
%
% Properties of OBJ_new now have values of 123 (for "prop1") and "hallo" (for "prop2")
%
% The class additionally implements a property "version", which I generally find very useful. It can be used while loading
% objects from disk to detect, if these are not up to date anymore. Since this property hidden and defineds as private,
% it may be redefined in any child-class and should not cause any problems.
%
% NOTE
% ====
% - This implementation (hopefully) also handles inherited properties, which are defined as 
%   "private" or "protected" in its superclass.
% - If properties contain child objects of the HANDLEPLUS-class, they are also deep copied.
%   If this is an unwanted behaviour, this needs to be changed manually.
% - This implementation DOES NOT handle cyclic references as discussed in this post:
%   <a href="matlab:web('mathforum.org/kb/message.jspa?messageID=7629086&tstart=0')"></a>.
%

% History
% =======
% 25.01.2012    Version 1 published to FEX
%
% Author
% ======
% Sebastian Hölz
% shoelz(at)geomar(dot)de
%    
    
    properties (SetAccess=private,GetAccess=private,Hidden)
        version=1;
    end
    
    % === Additional methods
    methods % (Sealed)  ... sealing these methods might make sense ...
        function obj_out = clone(obj, varargin)
            %
            % obj_new = OBJ.clone(varargin);
            %
            % Use this syntax to make a deep copy of an object OBJ, i.e. OBJ_OUT has the same field values, but will not behave as
            % a handle-copy of OBJ anymore. Additionally, varargin may contain valid property-value pairs, which will
            % change the copied values for the according property to the desired value.
            %
            obj_out = clone@handleplus(obj);
            if nargin>1; set(obj_out,varargin{:}); end            
            
        end
        function obj_out = copy(obj,varargin)
            %
            % obj_new = obj.copy; % Same as clone-method.
            %
            % Use this syntax to make a deep copy of an object OBJ, i.e. OBJ_OUT has the same field values,
            % but will not behave as a handle-copy of OBJ anymore. Additionally, varargin may contain valid property-value pairs, which will
            % change the copied values for the according property to the desired value.
            %
            obj_out = obj.clone(varargin{:});
        end
    end
    
end

