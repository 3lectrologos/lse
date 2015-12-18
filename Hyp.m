%--------------------------------------------------------------------------
% Author      : Alkis Gotovos <alkisg@inf.ethz.ch>
% Description : Hyperparameter abstract superclass
%--------------------------------------------------------------------------
classdef Hyp < handle
  
  properties (Access = public)
    fun; % Function forms (fun.mean, fun.cov, fun.lik)
    val; % Hyperparameter values (val.mean, val.cov, val.lik)
  end
  
  methods
    function self = Hyp(varargin)
      if nargin > 0
        self.fun = varargin{1};
      else
        self.fun.mean = @meanConst;
        self.fun.cov = @covSEiso;%{@covMaterniso, 5};
        self.fun.lik = @likGauss;
      end
      if nargin > 1
        self.val = varargin{2};
      else
        self.val.mean = 0;
        self.val.cov = [0; 0];
        self.val.lik = log(0.1);
      end
    end
  end
  
  methods (Abstract)
    % Update hyperparameters according to given training points
    update(obj, trainx, trainy);
  end
  
end