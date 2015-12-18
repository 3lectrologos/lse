%--------------------------------------------------------------------------
% Author      : Alkis Gotovos <alkisg@inf.ethz.ch>
% Description : Constant hyperparameters
%--------------------------------------------------------------------------
classdef Hyp_const < Hyp
  
  methods
    function self = Hyp_const(varargin)
      self = self@Hyp(varargin{:});
    end
    
    function update(self, trainx, trainy)
    end
  end
  
end