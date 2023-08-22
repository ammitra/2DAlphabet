'''make_pseudo.py

Module containing the functionality for producing a pseudo-data toy in the given regions of a 2DAlphabet fit,
using post-fit transfer functions to obtain the data-driven background estimates in the non-fail regions.

'''
import ROOT
import random

class PseudoData:
    '''Class to produce pseudo data from a given set of regions and postfit transfer functions.'''
    def __init__(self, twoD, regions=[], rpfs=[], findreplace={}):
	'''Creates a pseudodata toy in the given regions. Uses postfit transfer functions to obtain
	the data-driven background estimates in the non-fail regions, from which background PDFs are 
	generated and used to produce the toy data. 

	Args:
	    twoD (TwoDAlphabet): TwoDAlphabet object containing all information on the workspace
            regions (list(str)): List containing the names of the regions to be created, in the order
		in which the transfer functions would be applied. E.g., 'Fail', 'Loose', 'Pass'
            rpfs (list(TH2)): List containing the 2D histograms of the postfit transfer functions used, 
		in the order in which they were applied in the fit.
	    findreplace (dict): Non-nested dictionary with key-value pairs to find-replace in the TwoDAlphabet
		internal ledger.
	'''
	self.twoD = twoD
	self.ledger = self.twoD.ledger
	self.tag = twoD.tag
	self.regions = regions
	self.rpfs = rpfs
	# Create the dataframe with information on the nominal backgrounds and data files and histograms, 
	# either from the given regions or old regions if user needs to find and replace
	if findreplace:
	    # need to use the original regions (that will be replaced)
	    oldRegions = findreplace.keys()
	    self.df = self.select(oldRegions)
	    # Now find and replace the source histogram names from, e.g., control region to signal region
	    for find, replace in findreplace.items():
		#self.df['source_histname'] = self.df['source_histname'].apply(lambda x: x.replace(find, replace))
		self.df = self.df.replace(find, replace, regex=True) # Modifies both `region` and `source_histname`
	else:
	    # Check that a find-replace doesn't need to be performed
	    for i in range(len(self.df)):
		if not any([region in self.df.iloc[i].source_histname for region in self.regions]):
		    raise RuntimeError('Region names {} not found source in histogram {} for file {}. Run with the appropriate findreplace argument if you need to generate toys for regions which have not yet been fitted.'.format(self.regions, self.df.iloc[i].source_histname, self.df.iloc[i].source_filename))
	    # If all is good, create internal dataframe from the given regions list
	    self.df = self.df = self.select(self.regions)
	# Check that the number of regions and transfer functions makes sense
	if (len(regions)-len(rpfs)) != 1:
	    raise RuntimeError('There must be exactly 1 fewer transfer functions than regions between which to transfer.')
	# Get the output histogram names for the toy data in each region
	self.outHistNames = self._getOutputHistNames(self.df)
	self.nEvents = self._getNevents(self.df)

    def _select_nominal(self, row, args):
        '''Helper function to get only the Data and nominal background files+histos in the given fit regions'''
        regions = args[0]
        if row.region not in regions: return False
        elif row.process_type == 'DATA': return True
        elif row.process_type == 'BKG':
            if row.variation == 'nominal': return True
            else: return False
        else: return False

    def select(self, *args):
        '''Returns DataFrame containing the rows and columns selected by a lambda function'''
	def _select_nominal(row, args):
	    regions = args[0]
	    if row.region not in regions: return False
	    elif row.process_type == 'DATA': return True
	    elif row.process_type == 'BKG':
		if row.variation == 'nominal': return True
		else: return False
	    else: return False
        eval_lambda = lambda row: _select_nominal(row, args)
        df = self.ledger.df.loc[self.ledger.df.apply(eval_lambda,axis=1)]
        return df

    def _getOutputHistNames(self, df):
	'''Gets the histogram names from the find-replaced dataframe ledger for the toy data output histograms'''
	histnames = {} 
	for region, group in df.groupby('region'):
	    regionHistNames = group['source_histname'].drop_duplicates().to_list()
	    if len(regionHistNames) != 1:
		raise ValueError('There should only be one histogram name per region for data and nominal bkgs.')
	    histnames[region] = regionHistNames[0]
	return histnames

    def _getNevents(self, df):
        '''Searches through the internal dataframe and obtains the number of events IN DATA in each region.'''
	nEvents = {}
	for region, group in df.groupby('region'):
	    data = df.loc[df.region.eq(region) & df.process_type.eq('DATA')]
	    if (len(data) != 1):
		raise RuntimeError('There should only be one data histogram per region')
	    fName = data.iloc[0].source_filename
	    hName = data.iloc[0].source_histname
	    dataFile = ROOT.TFile.Open(fName,'READ')
	    dataHist = dataFile.Get(hName)
	    print('Obtaining total number of events in data in {}'.format(region))
	    n = dataHist.Integral()
	    nEvents[region] = n
	    dataFile.Close()
	return nEvents

    def GetFailTemplate(self):
	'''Subtract all MC background histograms from the data fail histogram and return the data-minus-bkg estimate.'''
	data_file = ROOT.TFile.Open(self.df.iloc[0].source_filename)
	data_hist = data_file.Get(self.df.iloc[0].source_histname)
	if (self.regions[0] not in data_hist.GetName()):
	    raise RuntimeError('Fail region {} not in data histogram name {}'.format(self.regions[0], data_hist.GetName()))
	bkg_est = data_hist.Clone('fail_template_bkg_est')
	fail_data = data_hist.Clone() # output the real fail data for future use
	fail_data.SetDirectory(0)
	bkg_est.SetDirectory(0)
	data_file.Close()
	print('Determining data minus background estimate in {}'.format(self.regions[0]))
	for region, group in self.df.groupby('region'):
	    if region != self.regions[0]: continue
	    bkg_sources = group.loc[group.process_type.eq('BKG')]
	    for i in range(len(bkg_sources)):
		bkg_file = ROOT.TFile.Open(bkg_sources.iloc[i].source_filename)
		bkg_hist = bkg_file.Get(bkg_sources.iloc[i].source_histname)
		print('Subtracting background {} from data in {}'.format(bkg_sources.iloc[i].process, region))
		bkg_est.Add(bkg_hist, -1)
		bkg_file.Close()
	return bkg_est, fail_data

    def GetBkgTemplates(self):
	'''Gets the summed backgrounds in each of the regions'''
	out = {}
        template_file = ROOT.TFile.Open(self.df.iloc[0].source_filename)
        template_hist = template_file.Get(self.df.iloc[0].source_histname)
	template_hist.SetDirectory(0)
	template_file.Close()
	for region in self.regions[1:]: # skip the fail region
	    print('Determining total MC background in {}'.format(region))
	    template = template_hist.Clone('TotalBkg_{}'.format(region))
	    template.Reset()
	    bkg_sources = self.df.loc[self.df.region.eq(region) & self.df.process_type.eq('BKG')]
	    for i in range(len(bkg_sources)):
                bkg_file = ROOT.TFile.Open(bkg_sources.iloc[i].source_filename)
                bkg_hist = bkg_file.Get(bkg_sources.iloc[i].source_histname)
		template.Add(bkg_hist)
		bkg_file.Close()
		out[region] = bkg_hist
	return out

    def GetAsimov(self, fail_template, bkg_templates):
	'''Constructs the data-driven background estimate in the non-fail regions from the
	product of the data-driven fail template and the transfer functions.

	Args:
	    fail_template (TH2): TH2 representing the data-bkg fail region estimate
	    bkg_templates (dict(TH2)): dictionary containing the key-value region name and 
		associated total MC background TH2

	Returns:
	    asimov (dict): key-value region name and associated total background (Asimov dataset)
	    transfers (dict): key-value region name and associated transfer function-based estimate
	'''
	# Find out how many transfers need to be performed
	nTransfers = len(self.regions)-1
	if nTransfers<1:
	    raise RuntimeError('There must be at least one transfer') # this should already be caught by now, but just in case
	# Get the total MC bkg distributions to which the transfer function-based estimates will be added
	asimov = bkg_templates.copy()
	transfers = {}

	for i, region in enumerate(self.regions[1:]): # skip the fail region
	    print('Performing transfer {}:\n\t{} -> {}'.format(i+1, self.regions[i], self.regions[i+1]))
            CurrentRegionIndex = self.regions.index(region)
            PreviousRegionIndex = CurrentRegionIndex - 1
	    PreviousRegionName = self.regions[PreviousRegionIndex]
	    TransferFuncName = self.rpfs[i].GetName()
	    outName = '{}_{}'.format(PreviousRegionName, TransferFuncName)
	    # Get the transfer-based estimate
	    if (i==0): # Performing the first transfer
		transferN = self.Multiply(fail_template, self.rpfs[i], outName)
	    else:
		CurrentRegionIndex = self.regions.index(region)
		PreviousRegionIndex = CurrentRegionIndex - 1
		transferN = self.Multiply(transfers[PreviousRegionName], self.rpfs[i], outName)
	    transferN.SetDirectory(0)
	    # Store the transfer-based estimate in case there is more than one transfer
	    transfers[region] = transferN
	    # Add the transfer-based bkg estimate to the total MC-based bkg estimate
	    asimov[region].Add(transferN)	

	return asimov, transfers

    def CreatePDF(self, asimovDatasets):
	'''Normalizes Asimov background distribution to unity to create a background PDF

	Args:
	    asimovDatasets (dict): key-value region name and associated total background distribution

	Returns:
	    PDFs (dict): key-value region name and associated total background PDF and number of events
		in data, in that region.
	'''
	PDFs = {}
	for region, asimovDataset in asimovDatasets.items():
	    PDF = asimovDatasets[region].Clone('pdf_{}'.format(region))
	    PDF.SetDirectory(0)
	    # The PDF must be scaled by the integral of teh *Asimov* dataset, but this is *NOT* the number of
	    # events that should be generated, since that is not representative of the number of events in the
	    # actual dataset, which is what we are trying to model
	    nEventsAsimov = PDF.Integral(1, PDF.GetNbinsX()+1, 1, PDF.GetNbinsY()+1)
	    PDF.Scale(1./nEventsAsimov)
	    PDFs[region] = [PDF, self.nEvents[region]] # return the actual number 
	return PDFs

    def CreateCDF(self, PDFs):
	'''Creates the cumulative distribution function for each region, based on the background PDF'''
	CDFs = {}
	for region, PDF in PDFs.items():
	    PDF = PDF[0] # just want the PDF histogram, not nEvents
	    nx = PDF.GetNbinsX()
	    ny = PDF.GetNbinsY()
	    CDF = ROOT.TH1F('CDF_{}'.format(region),"",nx*ny,0,nx*ny)
	    CDF.SetDirectory(0)
	    cumulativeBin = 0
	    for i in range(1,nx+1):
		for j in range(1,ny+1):
		    cumulativeBin += 1
		    PDFval = PDF.GetBinContent(i,j) + CDF.GetBinContent(cumulativeBin-1)
		    CDF.SetBinContent(cumulativeBin, PDFval)
	    # Check that the CDF is bounded within [0, 1)
	    if (CDF.GetMaximum() > 1.0):
		raise ValueError('CDF maximum ({}) > 1.0 - something has gone wrong'.format(CDF.GetMaximum()))
	    CDFs[region] = CDF
	return CDFs

    def _FindPDFintersection(self, rand, CDF):
	'''Returns the global bin corresponding to the CDF intersection'''
	found = False
	for i in range(1, CDF.GetNbinsX()+1):
	    PDFval = CDF.GetBinContent(i)
	    if (PDFval > rand):
		found = True
		return i
	if not found:
	    raise ValueError('Intersection with CDF not found, something is wrong. Make sure that the random number is bounded by the maximum value found in the CDF')

    def _GlobalBinTo2D(self, PDF, globalBin):
	'''Converts the global bin (from the CDF intersection) to the 2D PDF bin
	globalBins start from 1
	globalBin 1 = (1,1), 2 = (1,2), and so on
	'''
	globalBin = globalBin - 1
	NX = PDF.GetNbinsX()
	NY = PDF.GetNbinsY()
	localX = int(globalBin)/int(NY) + 1
	localY = globalBin % NY + 1
	return localX, localY

    def GeneratePseudoData(self, PDF, nEvents, CDF, name):
	'''Generates toy data from the CDF of a given region'''
	toy = PDF.Clone(name)
	toy.Reset()
	toy.SetDirectory(0)
	print('Generating {} events for toy data hist {}'.format(int(nEvents), name))
	for i in range(int(nEvents)):
	    rand = random.uniform(0,CDF.GetMaximum()) # constrain the random number to be within the max of the CDF
	    globalBin = self._FindPDFintersection(rand, CDF)
	    nx, ny = self._GlobalBinTo2D(PDF, globalBin)
	    nx = int(nx)
	    ny = int(ny)
	    toy.SetBinContent(nx, ny, toy.GetBinContent(nx,ny)+1)
	return toy

    def Multiply(self, h1, h2, name=''):
	'''Modified version of TH2::Multiply() which additionally zeroes negative bins'''
	#print('Multiplying {} x {}'.format(h1.GetName(),h2.GetName()))
	nh1 = (h1.GetNbinsX(),h1.GetNbinsY())
	nh2 = (h2.GetNbinsX(),h2.GetNbinsY())
	assert(nh1==nh2)
	#print('{} binning : {}\n{} binning : {}'.format(h1,nh1,h2,nh2))
	finalHist = h1.Clone(name)
	# loop over bins in h1
	for i in range(0, h1.GetNbinsX()+1):
	    for j in range(0, h1.GetNbinsY()+1):
		h1Val = h1.GetBinContent(i,j)
		if h1Val < 0:     # don't allow negative yields
		    h1Val = 0
		xVal = h1.GetXaxis().GetBinCenter(i)
		yVal = h1.GetYaxis().GetBinCenter(j)
		ih2 = h2.GetXaxis().FindBin(xVal)
		jh2 = h2.GetYaxis().FindBin(yVal)
		h2Val = h2.GetBinContent(ih2,jh2)
		finalHist.SetBinContent(i,j,h1Val*h2Val)
	finalHist.SetDirectory(0)
	return finalHist

    
