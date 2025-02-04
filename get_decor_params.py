import numbers
import warnings
import numpy as np
import argparse

# from https://github.com/nsmith-/rhalphalib/blob/master/rhalphalib/util.py
ROOFIT_HELPERS_INSTALLED = False
def install_roofit_helpers():
    global ROOFIT_HELPERS_INSTALLED
    if ROOFIT_HELPERS_INSTALLED:
        return
    ROOFIT_HELPERS_INSTALLED = True
    import ROOT as _ROOT

    def _RooFitResult_nameArray(self):
        """
        Returns a numpy array of floating parameter names
        """
        return np.array([p.GetName() for p in self.floatParsFinal()])

    _ROOT.RooFitResult.nameArray = _RooFitResult_nameArray

    def _RooFitResult_valueArray(self):
        """
        Returns a numpy array of floating parameter values
        """
        return np.array([p.getVal() for p in self.floatParsFinal()])

    _ROOT.RooFitResult.valueArray = _RooFitResult_valueArray

    def _RooFitResult_covarianceArray(self):
        """
        Returns a numpy array of floating parameter covariances
        """
        param_cov = self.covarianceMatrix()
        param_cov = np.frombuffer(param_cov.GetMatrixArray(), dtype="d", count=param_cov.GetNoElements())
        param_cov = param_cov.reshape(int(np.sqrt(param_cov.size)), -1)
        return param_cov

    _ROOT.RooFitResult.covarianceArray = _RooFitResult_covarianceArray


# from https://github.com/nsmith-/rhalphalib/blob/master/rhalphalib/parameter.pyclass Parameter(object):
class Parameter(object):
    def __init__(self, name, value):
        self._name = name
        self._value = value
        self._hasPrior = False
        self._intermediate = False

    def __repr__(self):
        return "<%s (%s) instance at 0x%x>" % (
            self.__class__.__name__,
            self._name,
            id(self),
        )

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        self._name = name

    @property
    def value(self):
        return self._value

    @property
    def intermediate(self):
        """
        An intermediate parameter is one that should not be explicitly rendered.
        The formula will be expanded recursively until it depends only on non-intermediate value.
        Only DependentParameters can be intermediate, hence one can modify this flag for them.
        """
        return self._intermediate

    def hasPrior(self):
        """
        True if the prior is not flat
        """
        return self._hasPrior

    @property
    def combinePrior(self):
        """
        By default assume param has no prior and we are just informing combine about it
        """
        return "flatParam"

    def getDependents(self, rendering=False, deep=False):
        return {self}

    def formula(self):
        return "{" + self._name + "}"

    def renderRoofit(self, workspace):
        raise NotImplementedError

    def _binary_op(self, opinfo, other):
        opname, op, right = opinfo
        if isinstance(other, Parameter):
            if right:
                name = other.name + opname + self.name
                out = DependentParameter(name, "{0}%s{1}" % op, other, self)
            else:
                name = self.name + opname + other.name
                out = DependentParameter(name, "{0}%s{1}" % op, self, other)
            out.intermediate = True
            return out
        elif isinstance(other, numbers.Number):
            if right:
                name = type(other).__name__ + opname + self.name
                out = DependentParameter(name, "%r%s{0}" % (other, op), self)
            else:
                name = self.name + opname + type(other).__name__
                out = DependentParameter(name, "{0}%s%r" % (op, other), self)
            out.intermediate = True
            return out
        return NotImplemented

    def __radd__(self, other):
        return self._binary_op(("_add_", "+", True), other)

    def __rsub__(self, other):
        return self._binary_op(("_sub_", "-", True), other)

    def __rmul__(self, other):
        return self._binary_op(("_mul_", "*", True), other)

    def __rtruediv__(self, other):
        return self._binary_op(("_div_", "/", True), other)

    def __rpow__(self, other):
        return self._binary_op(("_pow_", "**", True), other)

    def __add__(self, other):
        return self._binary_op(("_add_", "+", False), other)

    def __sub__(self, other):
        return self._binary_op(("_sub_", "-", False), other)

    def __mul__(self, other):
        return self._binary_op(("_mul_", "*", False), other)

    def __truediv__(self, other):
        return self._binary_op(("_div_", "/", False), other)

    def __pow__(self, other):
        return self._binary_op(("_pow_", "**", False), other)

    def max(self, val):
        """Return maximum out of param value and given ``val``"""
        return DependentParameter("max(%s,%s)" % (self.name, val), "TMath::Max({0}, %s)" % val, self)


class IndependentParameter(Parameter):
    DefaultRange = (-10, 10)

    def __init__(self, name, value, lo=None, hi=None, constant=False):
        super(IndependentParameter, self).__init__(name, value)
        self._lo = lo if lo is not None else self.DefaultRange[0]
        self._hi = hi if hi is not None else self.DefaultRange[1]
        self._constant = constant

    @Parameter.value.setter
    def value(self, val):
        self._value = val

    @property
    def lo(self):
        return self._lo

    @lo.setter
    def lo(self, lo):
        self._lo = lo

    @property
    def hi(self):
        return self._hi

    @hi.setter
    def hi(self, hi):
        self._hi = hi

    @property
    def constant(self):
        return self._constant

    @constant.setter
    def constant(self, const):
        self._constant = const

    def renderRoofit(self, workspace):
        import ROOT

        install_roofit_helpers()
        if workspace.var(self._name) == None:  # noqa: E711
            var = ROOT.RooRealVar(self._name, self._name, self._value, self._lo, self._hi)
            var.setAttribute("Constant", self._constant)
            workspace.add(var)
        return workspace.var(self._name)


class NuisanceParameter(IndependentParameter):
    def __init__(self, name, combinePrior, value=0, lo=None, hi=None):
        """
        A nuisance parameter.
        name: name of parameter
        combinePrior: one of 'shape', 'shapeN', 'lnN', etc.

        Render the prior somewhere else?  Probably in Model because the prior needs
        to be added at the RooSimultaneus level (I think)
        Filtering the set of model parameters for these classes can collect needed priors.
        """
        super(NuisanceParameter, self).__init__(name, value, lo, hi)
        self._hasPrior = True
        if combinePrior not in {"shape", "shapeN", "shapeU", "lnN", "lnU", "gmM", "trG", "param"}:
            raise ValueError("Unrecognized combine prior %s" % combinePrior)
        self._prior = combinePrior

    @property
    def combinePrior(self):
        return self._prior


class DependentParameter(Parameter):
    def __init__(self, name, formula, *dependents):
        """
        Create a dependent parameter
            name: name of parameter
            formula: a python format-string using only indices, e.g.
                '{0} + sin({1})*{2}'
        """
        super(DependentParameter, self).__init__(name, np.nan)
        if not all(isinstance(d, Parameter) for d in dependents):
            raise ValueError
        # TODO: validate formula for allowed functions
        self._formula = formula
        self._dependents = dependents

    @property
    def value(self):
        return eval(self.formula().format(**{p.name: p.value for p in self.getDependents(deep=True)}))

    @Parameter.intermediate.setter
    def intermediate(self, val):
        self._intermediate = val

    def getDependents(self, rendering=False, deep=False):
        """
        Return a set of parameters that this parameter depends on, which will be rendered.
        By default, this means all non-intermediate dependent parameters, recursively descending and stopping at
        the first renderable parameter (i.e. either non-intermediate or an IndependentParameter)
        If this parameter itself is renderable, we return a set of just this parameter.
        If rendering=True, we pass through this parameter if it is renderable.
        If deep=True, descend all the way to the IndependentParameters
        """
        dependents = set()
        if deep:
            for p in self._dependents:
                if isinstance(p, DependentParameter):
                    dependents.update(p.getDependents(deep=True))
                else:
                    dependents.add(p)
            return dependents
        if not (self.intermediate or rendering):
            return {self}
        for p in self._dependents:
            if p.intermediate:
                dependents.update(p.getDependents())
            else:
                dependents.add(p)
        return dependents

    def formula(self, rendering=False):
        if not (self.intermediate or rendering):
            return "{" + self.name + "}"
        return "(" + self._formula.format(*(p.formula() for p in self._dependents)) + ")"

    def renderRoofit(self, workspace):
        import ROOT

        install_roofit_helpers()
        if workspace.function(self._name) == None:  # noqa: E711
            if self.intermediate:
                # This is a warning because we should make sure the name does not conflict as
                # intermediate parameter names are often autogenerated and might not be unique/appropriate
                warnings.warn("Rendering intermediate parameter: %r" % self, RuntimeWarning)
                self.intermediate = False
            rooVars = [v.renderRoofit(workspace) for v in self.getDependents(rendering=True)]
            # Originally just passed the named variables to RooFormulaVar but it seems the TFormula class
            # is more sensitive to variable names than is reasonable, so we reindex here
            formula = self.formula(rendering=True).format(**{var.GetName(): "@%d" % i for i, var in enumerate(rooVars)})
            var = ROOT.RooFormulaVar(self._name, self._name, formula, ROOT.RooArgList.fromiter(rooVars))
            workspace.add(var)
        return workspace.function(self._name)


# from https://github.com/nsmith-/rhalphalib/blob/master/rhalphalib/function.py
class DecorrelatedNuisanceVector(object):
    def __init__(self, prefix, param_in, param_cov):
        if not isinstance(param_in, np.ndarray):
            raise ValueError("Expecting param_in to be numpy array")
        if not isinstance(param_cov, np.ndarray):
            raise ValueError("Expecting param_cov to be numpy array")
        if not (len(param_in.shape) == 1 and len(param_cov.shape) == 2 and param_cov.shape[0] == param_in.shape[0] and param_cov.shape[1] == param_in.shape[0]):
            raise ValueError("param_in and param_cov have mismatched shapes")

        _, s, v = np.linalg.svd(param_cov)
        self._transform = np.sqrt(s)[:, None] * v
        self._parameters = np.array([NuisanceParameter(prefix + str(i+1), "param") for i in range(param_in.size)])
        self._correlated = np.full(self._parameters.shape, None)
        self._correlated_str = []
        for i in range(self._parameters.size):
            coef = self._transform[:, i]
            order = np.argsort(np.abs(coef))
            self._correlated[i] = np.sum(self._parameters[order] * coef[order]) + param_in[i]
            self._correlated_str.append('('+' + '.join([f"{coef[o]:.3g}*{self._parameters[o].name}" for o in order])+f' + {param_in[i]:.3g})')

    @classmethod
    def fromRooFitResult(cls, prefix, fitresult, param_names=None):
        install_roofit_helpers()
        names = [p.GetName() for p in fitresult.floatParsFinal()]
        means = fitresult.valueArray()
        cov = fitresult.covarianceArray()
        if param_names is not None:
            pidx = np.array([names.index(pname) for pname in param_names])
            means = means[pidx]
            cov = cov[np.ix_(pidx, pidx)]
        out = cls(prefix, means, cov)
        if param_names is not None:
            for p, name in zip(out.correlated_params, param_names):
                p.name = name
        return out

    @property
    def parameters(self):
        return self._parameters

    @property
    def correlated_params(self):
        return self._correlated

    @property
    def correlated_str(self):
        return '\n'.join(self._correlated_str)


# from hessian.py
def get_objs(args):
    import ROOT

    fws, wsname = args.workspace.split(':')
    fin = ROOT.TFile.Open(fws)
    w = fin.Get(wsname)

    ffit, fitname = args.fit.split(':')
    fin2 = ROOT.TFile.Open(ffit)
    fit = fin2.Get(fitname)

    return w, fit

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Get decorrelated parameter expressions from covariance matrix",
    )
    subparsers = parser.add_subparsers()

    p_hessian = subparsers.add_parser('hessian')
    p_hessian.add_argument("-w", "--workspace", metavar="ROOTFILE:WORKSPACE", help="Workspace to load", required=True)
    p_hessian.add_argument("-f", "--fit", metavar="ROOTFILE:FIT_NAME", help="Fit result to load", required=True)
    p_hessian.add_argument("-m", "--model", help="Model to load", default="ModelConfig")
    p_hessian.set_defaults(mode='hessian')

    p_cov = subparsers.add_parser('cov')
    p_cov.add_argument("-t", "--tree", metavar="ROOTFILE:TREE_NAME", help="Tree to load", required=True)
    p_cov.add_argument("-p", "--parameters", nargs='+', help="Parameters to use", required=True)
    p_cov.set_defaults(mode='cov')

    args = parser.parse_args()

    if args.mode=='hessian':
        w, fit = get_objs(args)

        param_names = [p.GetName() for p in fit.floatParsFinal() if p.GetName().startswith("bsvj")]
        decoVector = DecorrelatedNuisanceVector.fromRooFitResult("@", fit, param_names)
    elif args.mode=='cov':
        import uproot as up
        fname, tname = args.tree.split(':')
        f = up.open(fname)
        pnames = args.parameters
        params = f[tname].arrays(pnames)

        # compute covariance matrix from scan in limit tree
        X = np.asarray([params[p] for p in pnames])
        means = np.mean(X, axis=1)
        cov = np.cov(X)
        decoVector = DecorrelatedNuisanceVector("@", means, cov)

    # output
    print(decoVector.correlated_str)

