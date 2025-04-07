"""
Scripts using building blocks in boosted_fits.py to create datacards
"""

import argparse, inspect, os, os.path as osp, re, json, itertools, sys, shutil
from pprint import pprint
from time import strftime, sleep
from copy import copy, deepcopy
from glob import glob

import numpy as np

import boosted_fits as bsvj

import ROOT # type: ignore
ROOT.RooMsgService.instance().setSilentMode(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.ERROR)

scripter = bsvj.Scripter()


@scripter
def plot_scipy_fits():
    parser = argparse.ArgumentParser(inspect.stack()[0][3])
    parser.add_argument('rootfile', type=str)
    parser.add_argument('-o', '--plotdir', type=str, default='plots_bkgfits_%b%d')
    parser.add_argument('-b', '--bdtcut', type=float, default=None)
    parser.add_argument('-n', '--npars', type=int, nargs='*')
    parser.add_argument('-p', '--pdftype', type=str, default=None, choices=bsvj.known_pdfs())

    args = parser.parse_args()
    plotdir = strftime(args.plotdir)
    if not osp.isdir(plotdir): os.makedirs(plotdir)

    import matplotlib.pyplot as plt # type: ignore
    bsvj.mpl_fontsizes()

    with bsvj.open_root(args.rootfile) as tf:

        def do_plot(tdir_name):
            tdir = tf.Get(tdir_name)
            bkg_hist = tdir.Get('Bkg')

            for pdf_type in bsvj.known_pdfs():
                if args.pdftype and pdf_type != args.pdftype: continue
                bsvj.logger.info('Fitting pdf_type=%s, tdir_name=%s', pdf_type, tdir_name)
                fig = plt.figure(figsize=(8,8))
                ax = fig.gca()
                binning, counts = bsvj.th1_binning_and_values(tdir.Get('Bkg'))
                bin_centers = np.array([.5*(l+r) for l, r in zip(binning[:-1], binning[1:])])
                # Bkg histogram
                ax.step(binning[:-1], counts, where='post', label='bkg {}'.format(tdir_name))
                # Fits
                if args.npars is not None and len(args.npars):
                    npars_iter = list(args.npars)
                else:
                    npars_iter = list(range(1,5) if pdf_type == 'alt' else range(2,6))
                for npars in npars_iter:
                    bsvj.logger.info('Fitting pdf_type=%s, tdir_name=%s, npars=%s', pdf_type, tdir_name, npars)
                    res = bsvj.fit_scipy(pdf_type, npars, bkg_hist)
                    y_pdf = bsvj.eval_expression(bsvj.pdf_expression(pdf_type, npars), [bin_centers] + list(res.x))
                    y_pdf = y_pdf/y_pdf.sum() * counts.sum()
                    chi2 = ((y_pdf-counts)**2 / y_pdf).sum() / (len(bin_centers) - npars)
                    label = '{}_npars{}, chi2={:.3f}, {}'.format(
                        pdf_type, npars, chi2,
                        ', '.join(['p{}={:.3f}'.format(i, v) for i, v in enumerate(res.x)])
                        )
                    ax.plot(bin_centers, y_pdf, label=label)
                ax.legend()
                ax.set_xlabel(r'$m_{T}$ (GeV)')
                ax.set_ylabel(r'$N_{events}$')
                ax.set_yscale('log')
                plt.savefig(osp.join(plotdir, tdir_name + '_' + pdf_type + '.png'), bbox_inches='tight')

        bdtcut = None
        if args.bdtcut is not None:
            tdir_name = 'bsvj_{:.1f}'.format(args.bdtcut).replace('.', 'p')
            do_plot(tdir_name)
        else:
            for tdir_name in [k.GetName() for k in tf.GetListOfKeys()]:
                do_plot(tdir_name)


@scripter
def plot_roofit_fits():
    parser = argparse.ArgumentParser(inspect.stack()[0][3])
    parser.add_argument('rootfile', type=str)
    parser.add_argument('-o', '--plotdir', type=str, default='plots_bkgfits_%b%d')
    parser.add_argument('-b', '--bdtcut', type=float, default=None)
    parser.add_argument('-n', '--npars', type=int, nargs='*')
    parser.add_argument('-p', '--pdftype', type=str, default=None, choices=bsvj.known_pdfs())
    args = parser.parse_args()
    plotdir = strftime(args.plotdir)
    if not osp.isdir(plotdir): os.makedirs(plotdir)

    with bsvj.open_root(args.rootfile) as tf:

        def do_plot(tdir_name):
            tdir = tf.Get(tdir_name)
            bkg_hist = tdir.Get('Bkg')
            mt = bsvj.get_mt_from_th1(bkg_hist)
            for pdf_type in bsvj.known_pdfs():
                if args.pdftype and pdf_type != args.pdftype: continue
                bsvj.logger.info('Fitting pdf_type=%s, tdir_name=%s', pdf_type, tdir_name)
                if args.npars is not None and len(args.npars):
                    npars_iter = list(args.npars)
                else:
                    npars_iter = list(range(1,5) if pdf_type == 'alt' else range(2,6))
                for npars in npars_iter:
                    bsvj.logger.info('Fitting pdf_type=%s, tdir_name=%s, npars=%s', pdf_type, tdir_name, npars)
                    res_scipy = bsvj.fit_scipy(pdf_type, npars, bkg_hist)
                    if len(res_scipy.x) != npars:
                        raise Exception(
                            'Wrong number of fitted parameters.'
                            ' Found {} parameters in scipy fit result, but npars is {}.'
                            ' Scipy fit result:\n{}'
                            .format(len(res_scipy.x), npars, res_scipy)
                            )
                    res_roofit_only = bsvj.fit_roofit(pdf_type, npars, bkg_hist)
                    res_roofit_wscipy = bsvj.fit_roofit(pdf_type, npars, bkg_hist, init_vals=res_scipy.x)
                    bsvj.plot_pdf_for_various_fitresults(
                        bsvj.make_pdf(pdf_type, npars, bkg_hist, mt=mt, name=bsvj.uid()),
                        [res_scipy, res_roofit_only, res_roofit_wscipy],
                        bsvj.th1_to_datahist(bkg_hist, mt=mt),
                        osp.join(plotdir, '{0}_{1}_npar{2}.png'.format(tdir.GetName(), pdf_type, npars)),
                        labels=['Scipy only', 'RooFit only', 'RooFit init w/ scipy']
                        )
                    print('-'*60)
                    print('Summary of varous fit strategies')
                    print('\nScipy only:')
                    print(res_scipy)
                    print('\nRooFit with initial parameters from Scipy:')
                    res_roofit_wscipy.Print()
                    print('\nRooFit with initial parameters set to 1.:')
                    res_roofit_only.Print()

        if args.bdtcut is not None:
            tdir_name = 'bsvj_{:.1f}'.format(args.bdtcut).replace('.', 'p')
            do_plot(tdir_name)
        else:
            for tdir_name in [k.GetName() for k in tf.GetListOfKeys()]:
                do_plot(tdir_name)


@scripter
def gen_datacards():
    jsons = bsvj.get_jsons()
    regions = bsvj.pull_arg('--regions', type=str, nargs='+').regions
    mtmin = bsvj.pull_arg('--mtmin', type=float, default=None).mtmin
    mtmax = bsvj.pull_arg('--mtmax', type=float, default=None).mtmax
    gof_type = bsvj.pull_arg('--gof-type', type=str, default='rss', choices=bsvj.choices('gof')).gof_type
    norm_type = bsvj.pull_arg('--norm-type', type=str, default='theta', choices=bsvj.choices('norm')).norm_type
    if mtmin is not None: jsons["mt_min"] = mtmin
    if mtmax is not None: jsons["mt_max"] = mtmax
    nosyst = bsvj.pull_arg('--nosyst', default=False, action="store_true").nosyst
    asimov_bkg = bsvj.pull_arg('--asimov-bkg', default=False, action="store_true").asimov_bkg
    winner = bsvj.pull_arg('--winner', default=None, nargs=2, action="append").winner
    winners = {a:int(b) for a,b in winner} if winner is not None else {}
    brute = bsvj.pull_arg('--brute', default=False, action="store_true").brute
    npar = bsvj.pull_arg('--npar', type=int, default=2).npar
    tf_from_mc = bsvj.pull_arg('--tf-from-mc', default=False, action="store_true").tf_from_mc
    suff = bsvj.pull_arg('--suff', type=str, default='').suff
    ftest = bsvj.pull_arg('--ftest', default=False, action="store_true").ftest
    ntoys = bsvj.read_arg('-t', type=int, default=0).t

    # use npar arg as max value in saturated gof f-test
    if ftest:
        if norm_type!="rhalpha":
            raise ValueError("F-test using saturated GoF only implemented for rhalpha method")
        npar_vals = range(npar)
    else:
        npar_vals = [npar]

    results = []
    for npar in npar_vals:
        npar_suff = suff
        if ftest:
            npar_name = f'npar{npar}'
            npar_suff = bsvj.joiner([suff,npar_name])
        input = bsvj.InputData(regions, norm_type, **jsons, asimov=asimov_bkg)
        dcfile = input.gen_datacard(nosyst=nosyst, gof_type=gof_type, winners=winners, brute=brute, npar=npar, tf_from_mc=tf_from_mc, suff=npar_suff)
        if ftest:
            result = {}
            # first use toys
            if ntoys>0:
                gof_file = gof_datacard(dcfile, npar_name, ntoys)
                result['toys'] = collect_gof(gof_file)
            # now use data
            with bsvj.reset_sys_argv():
                bsvj.pull_arg('-t', type=int)
                gof_data = gof_datacard(dcfile, npar_name)
                result['data'] = collect_gof(gof_data)
                if ntoys==0: result = result['data'][0]
            results.append((npar+1, result))

    # conduct ftest
    # todo: allow reusing existing result (--cache arg? or --npar -1?)
    if ftest:
        i_winner = bsvj.do_fisher_test(results, input.n_bins, a_crit=0.05, toys=ntoys>0)
        # assign i_winner as the "main" datacard
        outdir = osp.dirname(dcfile)
        npar_suff = f'_npar{i_winner}'
        wfiles = glob(osp.join(outdir,f'*{npar_suff}.*'))
        for wfile in wfiles:
            # leave intermediate f-test files for backup folder
            if osp.basename(wfile).startswith('higgsCombine'): continue
            shutil.move(wfile, wfile.replace(npar_suff,''))
        # put others into a backup folder
        bakdir = osp.join(outdir,'ftest')
        if not osp.isdir(bakdir): os.makedirs(bakdir)
        bfiles = glob(osp.join(outdir,'*_npar*.*'))
        for bfile in bfiles:
            shutil.move(bfile, osp.join(bakdir,osp.basename(bfile)))
        # also save results in backup folder for histogramming etc.
        with open(osp.join(bakdir,'results.py'), 'w') as resfile:
            resfile.write(f"winner = {i_winner}\n")
            resfile.write("results = "+repr(results)+"\n")


# to compute GoF statistic on datacard
def gof_datacard(dcfile, name, ntoys=0):
    # setup base command
    outdir = osp.dirname(dcfile)
    dc = bsvj.Datacard.from_txt(dcfile)
    cmd = bsvj.CombineCommand(dc)
    cmd.configure_from_command_line()
    cmd.name += osp.basename(dc.filename).replace('.txt','')

    # following sequence based on https://github.com/andrzejnovak/higgstocharm/blob/master/custom.py

    if ntoys!=0:
        # run GenerateOnly
        cmd_gen = bsvj.gen_toys(cmd)
        bsvj.run_combine_command(cmd_gen, outdir=outdir)
        cmd.kwargs['--toysFile'] = cmd_gen.outfile

    # run bestfit
    cmd_fit = bsvj.bestfit(cmd)
    cmd_fit.kwargs['--cminDefaultMinimizerStrategy'] = 0
    cmd_fit.kwargs['--robustFit'] = 1
    bsvj.run_combine_command(cmd_fit, outdir=outdir)

    # run saturated gof
    cmd_gof = cmd.copy()
    cmd_gof.method = 'GoodnessOfFit'
    cmd_gof.kwargs['--algo'] = 'saturated'
    cmd_gof.kwargs['--cminDefaultMinimizerStrategy'] = 0
    cmd_gof.kwargs['--snapshotName'] = 'MultiDimFit'
    cmd_gof.dc.filename = cmd_fit.outfile
    bsvj.run_combine_command(cmd_gof, outdir=outdir)

    return cmd_gof.outfile


# to process outputs from above for fisher test
def collect_gof(gof_file):
    # following https://github.com/andrzejnovak/higgstocharm/blob/master/new_plot_ftests.py
    out = {}
    try:
        file = ROOT.TFile.Open(gof_file)
        tree = file.Get("limit")
        for j in range(tree.GetEntries()):
            tree.GetEntry(j)
            # convert to nll
            out[j] = -2.0 * np.log(tree.limit)
    except:
        bsvj.logger.warning(f"Skipping {gof_file}")
        pass

    return out


def make_bestfit_and_scan_commands(txtfile, args=None):
    # this is shared between both commands, so pull it out first
    range = bsvj.pull_arg('-r', '--range', type=float, default=[-3., 5.], nargs=2).range
    if args is None: args = sys.argv[1:]
    with bsvj.set_args(sys.argv[:1] + args):
        dc = bsvj.Datacard.from_txt(txtfile)
        cmd = bsvj.CombineCommand(dc)
        cmd.name += osp.basename(dc.filename).replace('.txt','')
        scan = bsvj.scan(cmd, range)
        scan.name += 'Scan'
        scan.configure_from_command_line(scan=True)
        bestfit = bsvj.bestfit(cmd, range)
        bestfit.name += 'Bestfit'
        bestfit.configure_from_command_line()
    return bestfit, scan


@scripter
def bestfit(txtfile=None):
    range = bsvj.pull_arg('-r', '--range', type=float, default=None, nargs=2).range
    if txtfile is None:
        # Allow multiprocessing if multiple datacards are passed on the command line
        txtfiles = bsvj.pull_arg('datacards', type=str, nargs='+').datacards
        if len(txtfiles) > 1:
            # Call this function in a pool instead
            import multiprocessing as mp
            with mp.Manager() as manager:
                pool = mp.Pool(8)
                pool.map(bestfit, txtfiles)
                pool.close()
            return
        else:
            txtfile = osp.abspath(txtfiles[0])

    outdir = bsvj.pull_arg('-o', '--outdir', type=str, default=strftime('bestfits_%Y%m%d')).outdir
    outdir = osp.abspath(outdir)

    dc = bsvj.Datacard.from_txt(txtfile)
    cmd = bsvj.CombineCommand(dc)
    cmd.configure_from_command_line()
    cmd = bsvj.bestfit(cmd, range)
    cmd.raw = ' '.join(sys.argv[1:])
    cmd.name += 'Bestfit_' + osp.basename(txtfile).replace('.txt','')
    bsvj.run_combine_command(cmd, logfile=cmd.logfile, outdir=outdir)


@scripter
def gentoys():
    """
    Generate toys from datacards
    """
    datacards = bsvj.pull_arg('datacards', type=str, nargs='+').datacards
    outdir = bsvj.pull_arg('-o', '--outdir', type=str, default=strftime('toys_%Y%m%d')).outdir
    bsvj.logger.info(f'Output will be moved to {outdir}')

    for dc_file in datacards:
        dc = bsvj.Datacard.from_txt(dc_file)
        cmd = bsvj.CombineCommand(dc)
        cmd.configure_from_command_line()
        cmd.name += osp.basename(dc.filename).replace('.txt','')

        # Some specific settings for toy generation
        cmd = bsvj.gen_toys(cmd)

        bsvj.run_combine_command(cmd, outdir=outdir)


@scripter
def fittoys2():
    infiles = bsvj.pull_arg('infiles', type=str, nargs='+').infiles
    outdir = bsvj.pull_arg('-o', '--outdir', type=str, default=strftime('fittoys_%Y%m%d')).outdir
    bsvj.logger.info(f'Output will be moved to {outdir}')

    # Sort datacards and toysfiles
    datacards = []
    toysfiles = []
    for f in infiles:
        if 'GenerateOnly' in f:
            toysfiles.append(f)
        else:
            datacards.append(f)

    # Submit fit per datacard
    for dc_file in datacards:
        name = osp.basename(dc_file).replace('.txt','')
        dc = bsvj.Datacard.from_txt(dc_file)
        cmd = bsvj.CombineCommand(dc)
        cmd.configure_from_command_line()
        cmd = bsvj.bestfit(cmd)
        cmd.name += name

        for tf in toysfiles:
            if name in tf:
                cmd.kwargs['--toysFile'] = tf
                break
        else:
            raise Exception(
                f'Could not find a toy file for datacard {dc_file}; available toy files:\n'
                + "\n".join(toysfiles)
                )

        bsvj.run_combine_command(cmd, outdir=outdir)


@scripter
def fittoys():
    # cmdFit="combine ${DC_NAME_ALL}
    #    -M FitDiagnostics
    #    -n ${fitName}
    #    --toysFile higgsCombine${genName}.GenerateOnly.mH120.123456.root
    #    -t ${nTOYS}
    #    -v
    #    -1
    #    --toysFrequentist
    #    --saveToys
    #    --expectSignal ${expSig}
    #    --rMin ${rMin}
    #    --rMax ${rMax}
    #    --savePredictionsPerToy
    #    --bypassFrequentistFit
    #    --X-rtd MINIMIZER_MaxCalls=100000
    #    --setParameters $SetArgFitAll
    #    --freezeParameters $FrzArgFitAll
    #    --trackParameters $TrkArgFitAll"

    datacards = bsvj.pull_arg('datacards', type=str, nargs='+').datacards
    outdir = bsvj.pull_arg('-o', '--outdir', type=str, default=strftime('toyfits_%b%d')).outdir

    for dc_file in datacards:
        dc = bsvj.Datacard.from_txt(dc_file)
        cmd = bsvj.CombineCommand(dc)
        cmd.configure_from_command_line()
        cmd.name += osp.basename(dc.filename).replace('.txt','')

        cmd.method = 'FitDiagnostics'
        cmd.kwargs.pop('--algo', None)
        cmd.args.add('--toysFrequentist')
        cmd.args.add('--saveToys')
        cmd.args.add('--savePredictionsPerToy')
        cmd.args.add('--bypassFrequentistFit')
        cmd.kwargs['--X-rtd'] = 'MINIMIZER_MaxCalls=100000'

        toysFile = bsvj.pull_arg('--toysFile', required=True, type=str).toysFile

        if not '-t' in cmd.kwargs:
            with bsvj.open_root(toysFile) as f:
                cmd.kwargs['-t'] = f.Get('limit').GetEntries()

        assert '-t' in cmd.kwargs
        assert '--expectSignal' in cmd.kwargs

        bsvj.run_combine_command(cmd, outdir=outdir)

        fit_diag_file = 'fitDiagnostics{}.root'.format(cmd.name)
        os.rename(fit_diag_file, osp.join(outdir, fit_diag_file))

@scripter
def fithessian():
    datacards = bsvj.pull_arg('datacards', type=str, nargs='+').datacards
    outdir = bsvj.pull_arg('-o', '--outdir', type=str, default=strftime('hessianfits_%b%d')).outdir
    range = bsvj.pull_arg('-r', '--range', type=float, default=None, nargs=2).range

    for dc_file in datacards:
        dc = bsvj.Datacard.from_txt(dc_file)
        cmd = bsvj.CombineCommand(dc)
        cmd.configure_from_command_line()
        cmd.name += osp.basename(dc.filename).replace('.txt','')

        cmd.method = 'FitDiagnostics'
        cmd.kwargs.pop('--algo', None)
        cmd.args.add('--saveWorkspace')
        cmd.args.add('--saveShapes')
        cmd.kwargs['--X-rtd'] = ['REMOVE_CONSTANT_ZERO_POINT=1']
        cmd.args.add('--bypassFrequentistFit')
        cmd.kwargs['--X-rtd'].append('MINIMIZER_MaxCalls=100000')
        if range: cmd.add_range('r', range[0], range[1])

        assert '--expectSignal' in cmd.kwargs

        bsvj.run_combine_command(cmd, outdir=outdir)

        fit_diag_file = 'fitDiagnostics{}.root'.format(cmd.name)
        os.rename(fit_diag_file, osp.join(outdir, fit_diag_file))

@scripter
def impacts():
    dc_file = bsvj.pull_arg('datacard', type=str).datacard
    nfits = bsvj.pull_arg('--nfits', type=int, default=1).nfits

    dc = bsvj.Datacard.from_txt(osp.abspath(dc_file))
    base_cmd = bsvj.CombineCommand(dc)
    base_cmd.configure_from_command_line()
    if '--toysFrequentist' in base_cmd.args:
        base_cmd.args.remove('--toysFrequentist')

    workdir = strftime(f'impacts_cli_%Y%m%d_{osp.basename(dc_file).replace(".txt","")}')
    bsvj.logger.info(f'Executing from {workdir}')
    if not bsvj.DRYMODE:
        os.makedirs(workdir, exist_ok=True)
        os.chdir(workdir)

    # Initial fit
    cmd = base_cmd.copy()
    cmd.kwargs['--algo'] = 'singles'
    cmd.kwargs['--redefineSignalPOIs'] = 'r'
    cmd.kwargs['--floatOtherPOIs'] = 1
    cmd.kwargs['--saveInactivePOI'] = 1
    cmd.args.add('--saveWorkspace')
    cmd.name = '_initialFit_Test'
    if osp.isfile(cmd.outfile):
        bsvj.logger.warning(
            f'Initial fit output already exists, not running initial fit command: {cmd}'
            )
    else:
        bsvj.run_combine_command(cmd, logfile=cmd.logfile)
    initial_fit_outfile = cmd.outfile

    systs = []
    for syst in dc.syst_names_independent:
        if 'mcstat' in syst: continue
        if syst in base_cmd.freeze_parameters: continue
        systs.append(syst)
    bsvj.logger.info(f'Doing systematics: {" ".join(systs)}')

    # calculate all impacts
    combinetool_dofit_cmd = base_cmd.copy()
    combinetool_dofit_cmd.exe = 'combineTool.py'
    combinetool_dofit_cmd.method = 'Impacts'
    combinetool_dofit_cmd.dc.filename = osp.abspath(initial_fit_outfile)
    combinetool_dofit_cmd.kwargs['-m'] = 120
    combinetool_dofit_cmd.args.add('--doFits')
    combinetool_dofit_cmd.kwargs['--parallel'] = nfits
    combinetool_dofit_cmd.named.update(systs)
    combinetool_dofit_cmd.kwargs['--redefineSignalPOIs'] = 'r'
    combinetool_dofit_cmd.name = 'Test'
    bsvj.run_combine_command(combinetool_dofit_cmd, logfile=cmd.logfile)

    # combine all impacts
    impact_file = 'impacts.json'
    combinetool_json_cmd = (
        f'combineTool.py -M Impacts'
        f' -d {initial_fit_outfile} -m 120 -o {impact_file}'
        f' --named {",".join(systs)} --redefineSignalPOIs r'
        )
    bsvj.run_command(combinetool_json_cmd)

    # default impact plot w/ all systs
    plot_impacts_cmd = (
        f'plotImpacts.py'
        f' -i {impact_file} --label-size 0.047'
        f' -o impacts'
        )
    bsvj.run_command(plot_impacts_cmd)

    # modify impact json to remove bkg systs
    with open(impact_file,'r') as ifile:
        impacts_orig = json.load(ifile)
    impacts_nobkg = deepcopy(impacts_orig)
    impacts_nobkg["params"] = []
    for param in impacts_orig["params"]:
        if 'bkg' in param["name"].lower():
            continue
        else:
            impacts_nobkg["params"].append(param)
    # dump in the format used by combineTool/Impacts
    impact_file_nobkg = 'impacts_nobkg.json'
    with open(impact_file_nobkg,'w') as ofile:
        json.dump(impacts_nobkg, ofile, sort_keys=True, indent=2, separators=(',', ': '))
    plot_impacts_nobkg_cmd = (
        f'plotImpacts.py'
        f' -i {impact_file_nobkg} --label-size 0.047'
        f' -o impacts_nobkg'
        )
    bsvj.run_command(plot_impacts_nobkg_cmd)

@scripter
def likelihood_scan(args=None):
    """
    Runs a likelihood scan on a datacard
    """
    if args is None: args = sys.argv
    with bsvj.set_args(args):

        print(sys.argv)

        outdir = bsvj.pull_arg('-o', '--outdir', type=str).outdir
        txtfile = bsvj.pull_arg('datacard', type=str).datacard
        bestfit, scan = make_bestfit_and_scan_commands(txtfile)

        for cmd in [bestfit, scan]:
            bsvj.run_combine_command(cmd, logfile=cmd.logfile, outdir=outdir)
            if outdir is None:
                bsvj.logger.error('No outdir specified')


@scripter
def likelihood_scan_mp():
    """
    Like likelihood_scan, but accepts multiple datacards. 
    """
    datacards = bsvj.pull_arg('datacards', type=str, nargs='+').datacards
    outdir = bsvj.pull_arg('-o', '--outdir', type=str, default=strftime('scans_%Y%m%d')).outdir

    # Copy sys.argv per job, setting first argument to the datacard
    args = sys.argv[:]
    args.insert(1, datacards[0])
    args.extend(['--outdir', outdir])
    jobs = []
    for txtfile in datacards:
        args[1] = txtfile
        jobs.append(args[:])

    import multiprocessing
    p = multiprocessing.Pool(16)
    p.map(likelihood_scan, jobs)
    p.close()
    p.join()
    bsvj.logger.info('Finished pool')


@scripter
def printws():
    """
    Prints a workspace contents
    """
    parser = argparse.ArgumentParser(inspect.stack()[0][3])
    parser.add_argument('rootfile', type=str)
    parser.add_argument('-w', '--workspace', type=str)
    args = parser.parse_args()
    with bsvj.open_root(args.rootfile) as f:
        ws = bsvj.get_ws(f, args.workspace)
    ws.Print()
    return ws


@scripter
def remove_fsr():
    dc_files = bsvj.pull_arg('dcfiles', type=str, nargs='+').dcfiles
    for dc_file in dc_files:
        dc = bsvj.Datacard.from_txt(dc_file)
        dc.systs = [s for s in dc.systs if s[0] != 'fsr']
        bsvj.logger.info(f'Overwriting {dc_file}')
        with open(dc_file, 'w') as f:
            f.write(bsvj.parse_dc(dc))


if __name__ == '__main__':
    bsvj.debug(bsvj.pull_arg('-d', '--debug', action='store_true').debug)
    bsvj.drymode(bsvj.pull_arg('--dry', action='store_true').dry)
    scripter.run()
