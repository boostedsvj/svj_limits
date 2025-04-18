import argparse
import ROOT as r

r.gROOT.ProcessLine(r"""class RooFormula2 : public TNamed, public RooPrintable {
public:
RooArgList _origList;
std::vector<bool> _isCategory;
std::unique_ptr<TFormula> _tFormula;
std::string _reconstructFormula(std::string internalRepr) const {
  TString internalReprT = internalRepr.c_str();
  for (unsigned int i = 0; i < _origList.size(); ++i) {
    const auto& var = _origList[i];
    std::stringstream regexStr;
    regexStr << "x\\[" << i << "\\]|@" << i;
    TPRegexp regex(regexStr.str().c_str());
    std::string replacement = var.GetName();
    regex.Substitute(internalReprT, replacement.c_str());
  }
  return internalReprT.Data();
}
std::vector<std::string> _getVars() const {
  std::vector<std::string> results;
  results.reserve(_origList.size());
  for (unsigned int i = 0; i < _origList.size(); ++i) {
    results.push_back(_origList[i].GetName());
  }
  return results;
}
static std::string reconstructFormula(RooFormula* formula) { return reinterpret_cast<RooFormula2*>(formula)->_reconstructFormula(formula->GetTitle()); }
static std::vector<std::string> getVars(RooFormula* formula) { return reinterpret_cast<RooFormula2*>(formula)->_getVars(); }
};""")

def remove_dup(lines):
    seen = set()
    seen_add = seen.add
    return [x for x in lines if not (x in seen or seen_add(x))]

def get_deps(workspace, var):
    output = []
    ws_var = workspace.function(var)
    lhs = f"{ws_var.GetName()} ="
    if ws_var.ClassName()=="RooFormulaVar":
        formula = ws_var.formula()
        output.append(f"{lhs} {r.RooFormula2.reconstructFormula(formula)}")
        for dep in r.RooFormula2.getVars(formula):
            output.extend(get_deps(workspace, str(dep)))
    else:
        output.append(f"{lhs} {ws_var.getVal()}")
    return output

def main(filename, workspace, snapshot, var, verbose, all):
    file = r.TFile.Open(filename)
    ws = file.Get(workspace)
    if snapshot: ws.loadSnapshot(snapshot)
    lines = get_deps(ws, var)
    lines.reverse()
    lines = remove_dup(lines)
    lines_str = '\n'.join(lines)
    if verbose: print(lines_str)
    exec(lines_str)
    if all:
        for line in lines:
            v = line.split(" = ")[0]
            print(f"{v} = {eval(v)}")
    else:
        print(f"Value in workspace: {var} = {ws.function(var).getVal()}")
        print(f"    Value computed: {var} = {eval(var)}")

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", type=str, required=True, help="input filename")
    parser.add_argument("-w", "--workspace", type=str, required=True, help="input workspace name")
    parser.add_argument("-s", "--snapshot", type=str, help="snapshot to load")
    parser.add_argument("-x", "--var", type=str, required=True, help="variable/function name")
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="enable verbose printouts")
    parser.add_argument("-a", "--all", default=False, action="store_true", help="evaluate all variables")
    args = parser.parse_args()

    main(args.file, args.workspace, args.snapshot, args.var, args.verbose, args.all)
