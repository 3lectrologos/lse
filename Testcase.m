%--------------------------------------------------------------------------
% Author      : Alkis Gotovos <alkisg@inf.ethz.ch>
% Description : Testcase wrapper
%--------------------------------------------------------------------------
classdef Testcase
  
  properties
    name;
    samples;
    h;
    hyp;
    fun;
  end
  
  methods
    function self = Testcase(varargin)
      if nargin > 0
        self.samples = varargin{1};
      end
      if nargin > 1
        self.h = varargin{2};
      end
      if nargin > 2
        self.hyp = varargin{3};
      end
      if nargin > 3
        self.name = varargin{4};
      end
    end
  end
  
end