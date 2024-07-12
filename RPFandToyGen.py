from TwoDAlphabet.twoDalphabet import MakeCard, TwoDAlphabet
from TwoDAlphabet.alphawrap import BinnedDistribution, ParametricFunction, SemiParametricFunction
import ROOT
from argparse import ArgumentParser

def _generate_constraints(nparams):
    out = {}
    for i in range(nparams):
        if i == 0:
            out[i] = {"MIN":-500,"MAX":500}
        else:
            out[i] = {"MIN":-500,"MAX":500}
    return out

_rpf_options = {
    '0x0': {
        'form': '(@0)',
        'constraints': _generate_constraints(1)
    },
    '1x0': {
        'form': '(@0+@1*x)',
        'constraints': _generate_constraints(2)
    },
    '0x1': {
        'form': '(@0+@1*y)',
        'constraints': _generate_constraints(2)
    },
    '1x1': {
        'form': '(@0+@1*x)*(@2+@3*y)',
        'constraints': _generate_constraints(3)
    },
    '2x1': {
        'form': '(@0+@1*x+@2*x**2)*(@3+@4*y)',
        'constraints': _generate_constraints(4)
    },
    '2x2': {
        'form': '(@0+@1*x+@2*x**2)*(@3+@4*y*@5*y**2)',
        'constraints': _generate_constraints(4)
    },
    '3x2': {
        'form': '(@0+@1*x+@2*x**2+@3*x**3)*(@4+@5*y)',
        'constraints': _generate_constraints(4)
    }
}

parser = ArgumentParser()
parser.add_argument('-w', type=str, dest='workspace',
                    action='store', required=True,
                    help='Setname to process.')
parser.add_argument('-f', type=str, dest='func',
                    action='store', required=True,
                    help='transfer function used. e.g. "b_2x1"')
args = parser.parse_args()

workspace = args.workspace

# Load up a previous workspace with fits already performed
twoD = TwoDAlphabet(workspace,'{}/runConfig.json'.format(workspace),findreplace={},loadPrevious=True)

# The new function saves the transfer function shapes (in the original input binning) to a
# new ROOT file, and returns the histograms themselves for use.
rpf_dict = twoD.GetTransferFunctionShapes(binName='default',rpfOpts=_rpf_options)

# Now we produce a pseudo-data toy in the SR.
# First, declare which regions we want to create a toy for. 
regions = ['SR_loose','SR_pass'] 
# Note that these regions may not be included in the TwoDAlphabet's configuration JSON, 
# but do exist in the input datasets as histograms. As such, we need to do a quick find-and-replace
findreplace = {'CR_loose':'SR_loose', 'CR_pass':'SR_pass'}
# We also need to tell the toy generation method which postfit transfer function to use to produce the QCD estimate 
# in all the non-fail regions. In this case, we are only fitting Loose -> Pass, so there should only be one function
rpfs = [rpf_dict[args.func]]
# Now we can pass this all to the toy generation method, which will output a root file containing the histograms for:
# 	- SR_loose (we want this to be blinded for now, since it's not entirely signal-free)
#	- SR_pass
twoD.MakePseudoData(regions, rpfs, findreplace, blindFail=True, poisson=True)

