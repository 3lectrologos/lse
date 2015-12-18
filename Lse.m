%--------------------------------------------------------------------------
% Author      : Alkis Gotovos <alkisg@inf.ethz.ch>
% Description : LSE implementation
%--------------------------------------------------------------------------
classdef Lse < handleplus
  
  properties
    h;        % Threshold
    epsilon;  % Accuracy parameter
    samples;  % Available sample points (samples.x, samples.y)
    gtruth;   % Ground truth
    ns;       % Number of available sample points (length(samples.y))
    fun;      % Function to compute ys
    vt;       % Visited set
    xvt;      % Extra visited points during traversal of planned path
    ht;       % Higher set
    lt;       % Lower set
    ut;       % Unknown set
    mt;       % Maximum set (subset of ut)
    om;       % Omega parameter for implicit threhold version
    rlo; rhi; % Current lower and upper bound of each conf. interval
    t;        % Current time step
    bval;     % If bval == -1, use method, else use given value
    reclass;  % Perform reclassification after termination?
    pshow;    % Show plot?
    mdl;      % GP model
    rule;     % Next measurement selection rule
    pedantic; % Strictly stick to algorithms as described in theory
  end
  
  methods
    function self = Lse(varargin)
      if nargin > 1
        testcase = varargin{1};
        self.samples = testcase.samples;
        if isfield(testcase.samples, 'gt')
          self.gtruth = testcase.samples.gt;
        elseif isfield(testcase.samples, 'y')
          self.gtruth = testcase.samples.y;
        end
        if isprop(testcase, 'fun')
          self.fun = testcase.fun;
        end
        self.h = testcase.h;
        self.mdl = Model(size(testcase.samples.x, 2), testcase.hyp);
        self.epsilon = varargin{2};
        self.rule = 'amb';
        self.reclass = true;
        self.pshow = true;
        self.bval = -1;
        self.pedantic = false;
        self.reset();
      end
    end
    
    function reset(self)
      n = size(self.samples.x, 1);
      self.ns = n;
      if ~isfield(self.samples, 'y') || isempty(self.samples.y)
        if ~isprop(self, 'fun')
          error('Need either samples.y or fun');
        else
          self.samples.y = NaN*ones(n, 1);
        end
      end
      self.ht = [];
      self.lt = [];
      self.ut = 1:n;
      self.vt = [];
      self.xvt = [];
      self.mt = [];
      self.rlo = -Inf*ones(n, 1);
      self.rhi = Inf*ones(n, 1);
      if ~isempty(self.om)
        self.h = [];
      end
      self.t = 1;
      self.mdl.reset();
    end
    
    function y = gety(self, idx)
      nanidx = isnan(self.samples.y(idx));
      if ~isempty(idx(nanidx))
        res = self.fun(self.samples.x(idx(nanidx), :));
        self.samples.y(idx(nanidx)) = res;
      end
      y = self.samples.y(idx);
    end
    
    function res = b(self)
      if self.bval == -1
        %delta = 0.01;
        %res = (2*log(self.ns*(pi^2)*(self.t^2)/(6*delta)))^(0.5);
        res = 3;
      else
        res = self.bval;
      end
    end
    
    function [val, idx] = rmax(~, x)
      val = max(x);
      idxs = find(x == max(x));
      idx = idxs(randsample(length(idxs), 1));
    end
    
    function nextidx = next(self)
      assert(~isempty(self.rlo));
      if strcmp(self.rule, 'var')
        [~, nextidx] = self.rmax(self.rhi - self.rlo);
      elseif strcmp(self.rule, 'amb')
        hm = self.geth();
        d1 = self.rhi - hm;
        d2 = hm - self.rlo;
        [~, nextidx] = self.rmax(min([d1 d2], [], 2));
      elseif strcmp(self.rule, 'random')
        nextidx = randsample(1:length(self.rlo), 1);
      else
        error('Unknown selection rule for next measurement.');
      end
    end
    
    function [idx, nextidx] = get_next(self)
      hm = self.geth();
      nextidx = self.next();
      idx = self.ut(nextidx);
      self.vt = [self.vt idx];
      self.ut(nextidx) = [];
      self.mt(self.mt == idx) = [];
      self.rlo(nextidx) = [];
      self.rhi(nextidx) = [];
      if isinf(self.gety(idx))
        [idx, nextidx] = self.get_next();
      elseif self.gety(idx) > hm
        self.ht = [self.ht idx];
      else
        self.lt = [self.lt idx];
      end
    end
    
    function [a, b, c, d] = geth(self)
      if ~isempty(self.om)
        [~, ~, m, v] = self.mdl.inf(self.samples.x);
        a = self.om*max(m);
        if self.pedantic
          b = self.om*max(self.rlo);
          c = self.om*max(self.rhi);
        else
          b = self.om*max(m-self.b*sqrt(v));
          c = self.om*max(m+self.b*sqrt(v));
        end
        d = max(m-self.b*sqrt(v));
      else
        a = self.h;
        b = self.h;
        c = self.h;
        d = [];
      end
    end
    
    function classify(self)
      [~, hl, hh, ml] = self.geth();
      lo = self.rhi <= hl + self.epsilon;
      hi = self.rlo >= hh - self.epsilon;
      % In the implicit case, self.ut is what we call zt in the formal alg.
      if ~isempty(self.om) && ~isempty(self.ut)
        if self.pedantic
          tut = setdiff(self.ut, self.mt);
          nopt = self.rhi < ml;
          self.mt = union(self.mt,...
                          intersect(tut, self.ut(~nopt & (lo|hi))));
          [~, idxs] = intersect(self.ut, self.mt, 'stable');
          idxm = false(length(self.ut), 1);
          idxm(idxs) = true;
          hi = nopt & hi & (~idxm);
          lo = nopt & lo & (~idxm);
        else
          nopt = self.rhi < ml;
          self.mt = self.ut(~nopt & (lo|hi));
          hi = nopt & hi;
          lo = nopt & lo;
        end
      end
      self.lt = reshape(union(self.lt, self.ut(lo)), 1, []);
      self.ht = reshape(union(self.ht, self.ut(hi)), 1, []);
      self.ut(lo|hi) = [];
      self.rlo(lo|hi) = [];
      self.rhi(lo|hi) = [];
    end
    
    function strict_classify(self)
      hm = self.geth();
      [~, ~, m] = self.mdl.inf(self.samples.x);
      idx = 1:length(self.samples.y);
      self.ht = idx(m > hm);
      self.lt = idx(m <= hm);
      self.ut = [];
      self.mt = [];
      self.rhi = [];
      self.rlo = [];
    end
    
    function [m, v] = inf(self)
      [~, ~, m, v] = self.mdl.inf(self.samples.x(self.ut, :));
      if self.pedantic
        self.rlo = max([m - self.b*sqrt(v) self.rlo], [], 2);
        self.rhi = min([m + self.b*sqrt(v) self.rhi], [], 2);
      else
        self.rlo = m - self.b*sqrt(v);
        self.rhi = m + self.b*sqrt(v);
      end
    end
    
    function [m, v] = inf_no_is(self)
      [~, ~, m, v] = self.mdl.inf(self.samples.x(self.ut, :));
      self.rlo = m - self.b*sqrt(v);
      self.rhi = m + self.b*sqrt(v);
    end
    
    function run(self, varargin)
      if nargin > 1
        maxiter = varargin{1};
      else
        maxiter = Inf;
      end
      self.reset();
      % Initialization
      ninit = 0;
      if ninit > 0
        srule = self.rule;
        self.rule = 'var';
        for i = 1:min(ninit, maxiter)
          self.inf();
          self.plot; drawnow;
          idx = self.get_next();
          self.mdl.add(self.samples.x(idx, :), self.gety(idx));
          self.t = self.t + 1;
        end
        self.rule = srule;
      end
      while ~isempty(setdiff(self.ut, self.mt)) && self.t <= maxiter
        self.inf();
        self.classify();
        if isempty(setdiff(self.ut, self.mt))
          break;
        end
        idx = self.get_next();
        self.mdl.add(self.samples.x(idx, :), self.gety(idx));
        self.t = self.t + 1;
        self.plot; drawnow;
      end
      if self.reclass
        self.strict_classify();
        self.plot; drawnow;
      end
    end
    
    function run_batch(self, nbatch, varargin)
      assert(nbatch > 0);
      if nargin > 2 && ~isempty(varargin{1})
        ninit = varargin{1};
      end
      if nargin > 3 && ~isempty(varargin{2})
        maxiter = varargin{2};
      else
        maxiter = Inf;
      end
      self.reset();
      % Initialization
      srule = self.rule;
      self.rule = 'var';
      for i = 1:min(ninit, maxiter)
        self.inf();
        self.plot; drawnow;
        idx = self.get_next();
        self.mdl.add(self.samples.x(idx, :), self.gety(idx));
        self.t = self.t + 1;
      end
      % Batch
      self.rule = srule;
      self.classify();
      while ~isempty(setdiff(self.ut, self.mt)) && self.t <= maxiter
        for i = 1:nbatch
          m = self.inf_no_is();
          if isempty(setdiff(self.ut, self.mt))
            i = i - 1;
            break;
          end
          [idx, utidx] = self.get_next();
          self.mdl.add(self.samples.x(idx, :), m(utidx));
          if isempty(self.ut)
            break;
          end
          self.t = self.t + 1;
        end
        % Add actual measurements of batch-selected points to model
        batch = self.vt(end-i+1:end);
        self.mdl.remove(i);
        self.mdl.add(self.samples.x(batch, :), self.gety(batch));
        if isempty(setdiff(self.ut, self.mt))
          self.plot; drawnow;
          break;
        end
        self.inf();
        self.classify();
        self.plot; drawnow;
      end
      if self.reclass
        self.strict_classify();
        self.plot; drawnow;
      end
    end
    
    % Return F1-score and classification error
    function [fsc, cerr] = eval(self)
      if ~isempty(self.om)
        hinf = self.om*max(self.gtruth);
      else
        hinf = self.h;
      end
      self.classify();
      tp = self.ht(self.gtruth(self.ht) > hinf);
      fp = self.ht(self.gtruth(self.ht) < hinf);
      fn = self.lt(self.gtruth(self.lt) > hinf);
      ntp = length(tp);
      nfp = length(fp);
      nfn = length(fn);
      prec = ntp/(ntp + nfp + eps);
      rec = ntp/(ntp + nfn + eps);
      % F1-score
      fsc = 2*prec*rec/(prec + rec + eps);
      % Classification error
      cerr = sum(abs(self.gtruth([fp fn]) - hinf));
    end
    
    function plot(self)
      if ~self.pshow
        return;
      end
      hdl = newplot();
      flag = false;
      if ~ishold()
        flag = true;
        hold on;
      end
      dcol = [0.2 0.2 0.2];
      dlcol = [0.5 0.5 0.5];
      lcol = [1 1 1];
      bcol = [0 0.5 0.5];
      rcol = [1 0.5 0.3];
      mcol = [0 1 0.2];
      plot(hdl, self.samples.x(self.ht, 1), self.samples.x(self.ht, 2),...
           'o', 'MarkerFaceColor', rcol, 'MarkerEdgeColor', rcol,...
           'MarkerSize', 5);
      plot(hdl, self.samples.x(self.lt, 1), self.samples.x(self.lt, 2),...
           'o', 'MarkerFaceColor', bcol, 'MarkerEdgeColor', bcol,...
           'MarkerSize', 5);
      plot(hdl, self.samples.x(self.ut, 1), self.samples.x(self.ut, 2),...
           'o', 'MarkerFaceColor', dlcol, 'MarkerEdgeColor', dlcol,...
           'MarkerSize', 4);
      plot(hdl, self.samples.x(self.mt, 1), self.samples.x(self.mt, 2),...
           'o', 'MarkerFaceColor', mcol, 'MarkerEdgeColor', mcol,...
           'MarkerSize', 4);
      if ~isempty(self.xvt)
        plot(hdl, self.xvt(:, 1), self.xvt(:, 2), 'x', 'LineWidth', 2,...
             'MarkerFaceColor', dcol, 'MarkerEdgeColor', dcol,...
             'MarkerSize', 4);
      end
      plot(hdl, self.samples.x(self.vt, 1), self.samples.x(self.vt, 2),...
           'o', 'MarkerFaceColor', dcol, 'MarkerEdgeColor', dcol,...
           'MarkerSize', 8);
      vl = intersect(self.vt, self.lt);
      plot(hdl, self.samples.x(vl, 1), self.samples.x(vl, 2), 'o',...
           'MarkerFaceColor', bcol, 'MarkerEdgeColor', bcol,...
           'MarkerSize', 6);
      vh = intersect(self.vt, self.ht);
      plot(hdl, self.samples.x(vh, 1), self.samples.x(vh, 2), 'o',...
           'MarkerFaceColor', rcol, 'MarkerEdgeColor', rcol,...
           'MarkerSize', 6);
      plot(hdl, self.samples.x(self.vt, 1), self.samples.x(self.vt, 2),...
           'o', 'MarkerFaceColor', lcol, 'MarkerEdgeColor', lcol,...
           'MarkerSize', 3);
      title(sprintf('Samples = %d', length(self.vt)));
      if flag
        hold off;
      end
    end
  end
  
end
