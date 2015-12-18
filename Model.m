%--------------------------------------------------------------------------
% Author      : Alkis Gotovos <alkisg@inf.ethz.ch>
% Description : GP model class, providing funtionality for adding/removing
%               training samples, doing  GP inference, and fitting
%               hyperparameters to given training samples
%--------------------------------------------------------------------------
classdef Model < handleplus
  
  properties
    d;       % Input space dimension
    train;   % Measurements (train.x, train.y)
    hyp;     % Model hyperparameters
  end
  
  methods
    function self = Model(varargin)
      if nargin > 1
        self.d = varargin{1};
        self.hyp = varargin{2};
      end
      self.reset();
    end
    
    function reset(self)
      self.train.x = zeros(0, self.d);
      self.train.y = zeros(0, 1);
    end
    
    % Add one or more measurements
    function add(self, x, y)
      if isempty(x) && isempty(y)
        return;
      end
      assert(size(x, 1) == size(y, 1));
      assert(size(x, 2) == self.d);
      assert(size(y, 2) == 1);
      self.train.x = [self.train.x; x];
      self.train.y = [self.train.y; y];
    end
    
    % Remove the n most recent measurements
    function remove(self, n)
      assert(n >= 0);
      self.train.x(end-n+1:end, :) = [];
      self.train.y(end-n+1:end) = [];
    end
    
    % GP model inference for given test points
    function [my, vy, mf, vf] = inf(self, testx)
      [my, vy, mf, vf] = gp(self.hyp.val, @infExact, self.hyp.fun.mean,...
                            self.hyp.fun.cov, self.hyp.fun.lik,...
                            self.train.x, self.train.y, testx);
    end
    
    % Conditional mutual info between GP and testx given the existing
    % training points
    function mi = mutual_info(self, testx)
      if ~iscell(self.hyp.fun.cov)
        self.hyp.fun.cov = {self.hyp.fun.cov};
      end
      Kbb = feval(self.hyp.fun.cov{:}, self.hyp.val.cov, testx);
      Kab = feval(self.hyp.fun.cov{:}, self.hyp.val.cov, self.train.x,...
                  testx);
      Kaa = feval(self.hyp.fun.cov{:}, self.hyp.val.cov, self.train.x);
      sn2 = exp(2*self.hyp.val.lik);
      na = length(Kaa);
      S = Kbb - Kab'*((Kaa + sn2*eye(na))\Kab);
      ns = length(S);
      mi = 0.5*logdet(eye(ns) + sn2^(-1)*S);
    end
    
    % GP model inference on grid
    function [my, vy, mf, vf] = inf_grid(self, x1, x2)
      n = length(x1);
      x = [x1(:) x2(:)];
      [mry, vry, mrf, vrf] = self.inf(x);
      my = reshape(mry, n, n);
      vy = reshape(vry, n, n);
      mf = reshape(mrf, n, n);
      vf = reshape(vrf, n, n);
    end
    
    % Update model hyperparameters
    function update_hyp(self)
      self.hyp.update(self.train.x, self.train.y);
    end
  end
end