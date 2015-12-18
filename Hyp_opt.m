%--------------------------------------------------------------------------
% Author      : Alkis Gotovos <alkisg@inf.ethz.ch>
% Description : Optimized hyperparameters
%--------------------------------------------------------------------------
classdef Hyp_opt < Hyp
  
  methods
    function self = Hyp_opt(varargin)
      self = self@Hyp(varargin{:});
    end
    
    function update(obj, trainx, trainy)
      obj.val = minimize(obj.val, @gp, -100, @infExact,...
                         obj.fun.mean, obj.fun.cov, obj.fun.lik,...
                         trainx, trainy);
    end
  end
  
end