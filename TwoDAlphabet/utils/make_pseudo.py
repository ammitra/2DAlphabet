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
	# Create an empty template for filling purposes
	temp_file = ROOT.TFile.Open(self.df.iloc[0].source_filename)
        temp_hist = temp_file.Get(self.df.iloc[0].source_histname)
	template_hist = temp_hist.Clone('template_hist')
	template_hist.Reset()
	self.template_hist = template_hist
	self.template_hist.SetDirectory(0)
	temp_file.Close()

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

    def GetData(self, region):
	'''Get the data histogram for a given region
	Returns:
	    region_data (TH2): actual data distribution in given region
	'''
	region_data = None
	for r, group in self.df.groupby('region'):
	    if r != region: continue
	    print('Obtaining data in {}'.format(region))
	    data_source = group.loc[group.process_type.eq('DATA')]
	    data_file = ROOT.TFile.Open(data_source.iloc[0].source_filename)
	    data_hist = data_file.Get(data_source.iloc[0].source_histname)
	    region_data = data_hist.Clone('ACTUAL_DATA_{}'.format(region))
	region_data.SetDirectory(0) # will act as an exception if something broken
	return region_data

    def GetBackgrounds(self, region):
	'''Get a histogram of all nominal backgrounds summed for a given region
	Returns:
	    total_bkg (TH2): total background distribution in a given region
	'''
	print('Obtaining total background in {}'.format(region))
	# Make a blank template to store the total bkg for the given region
	total_bkg = self.template_hist.Clone('TotalMCbkg_{}'.format(region))
	# Obtain all the background files for the given region
	for r, group in self.df.groupby('region'):
	    if r != region: continue
	    bkg_sources = group.loc[group.process_type.eq('BKG')]
	    for i in range(len(bkg_sources)):
		bkg_file = ROOT.TFile.Open(bkg_sources.iloc[i].source_filename)
                bkg_hist = bkg_file.Get(bkg_sources.iloc[i].source_histname)
		print('Adding background {} to total bkg in {}'.format(bkg_sources.iloc[i].process, region))
		total_bkg.Add(bkg_hist)
		bkg_file.Close()
	total_bkg.SetDirectory(0)
	return total_bkg

    def GetDataMinusBackgrounds(self, region, data, nominalMC):
	'''Get the histogram of the data minus nominalMC backgrounds for a given region
	Args:
	    region (str): region to obtain data-nominalMC bkg for
	    data (TH2): data in the given region
	    nominalMC (TH2): nominal MC backgrounds in the given region
	Returns:
	    dataMinusBkg (TH2): data minus nominalMC background for a given region
	'''
	dataMinusBkg = self.template_hist.Clone('DataMinusMCBkg_{}'.format(region))
	dataMinusBkg.Add(data)
	dataMinusBkg.Add(nominalMC, -1.)
	return dataMinusBkg

    def GetDataDrivenPrediction(self, dataMinusBkgFail):
	'''Get the data-driven prediction for the non-fail regions using the transfer functions
	and the data-minus-nominalMC distribution in the fail.
	Args:
	    dataMinusBkgFail (TH2): The data-minus-nominalMC distribution in the fail, used to 
		begin the chain of transfer functions
	Returns:
	    dataDriven (dict(TH2)): key-value dict containing the region name and the data-driven
		prediction in that region
	'''
	# Find out how many transfers need to be performed
	nTransfers = len(self.regions)-1
	if nTransfers < 1:
	    raise RuntimeError('There must be at least one transfer') # this should already be caught by now, but just in case
	dataDriven = {}
	for i, region in enumerate(self.regions[1:]): # skip the fail region
	    print('Performing transfer {}:\n\t{} -> {}'.format(i+1, self.regions[i], self.regions[i+1]))
	    currentRegionIdx = self.regions.index(region)
	    previousRegionIdx = currentRegionIdx - 1
	    previousRegionName = self.regions[previousRegionIdx]
	    transferFuncName = self.rpfs[i].GetName()
	    outName = '{}_times_{}_equals_{}'.format(previousRegionName, transferFuncName, self.regions[currentRegionIdx])
	    # Get the transfer-based estimate
	    if (i==0):
		transferN = self.Multiply(dataMinusBkgFail, self.rpfs[i], outName)
	    else:
		transferN = self.Multiply(dataDriven[previousRegionName], self.rpfs[i], outName)
	    transferN.SetDirectory(0)
	    # Store the transfer-based estimate in case there is more than one transfer
	    dataDriven[region] = transferN
	# Return the output dict
	return dataDriven

    def GetAsimov(self, region, nominalMC, dataDriven):
	'''Get the Asimov dataset (nominalMC + data-driven prediction) for given region
        Args:
            region (str): name of region to create PDF for
            nominalMC (TH2): summed nominal MC backgrounds for given region
            dataDriven (TH2): data-driven prediction for the given region. If this is the fail region, then the
                data-driven prediction should be given by the data-minus-nominalMC distribution in the fail.
                Otherwise, the data-driven prediction should be given by that which was obtained through successive
                applications of the transfer function. See GetDataDrivenPrediction()
        Returns:
            asimov (TH2): The asimov dataset for the given region
	'''
	asimov = self.template_hist.Clone('Asimov_{}'.format(region))
	asimov.Add(nominalMC)
	asimov.Add(dataDriven)
	asimov.SetDirectory(0)
	return asimov

    def GetPDF(self, region, asimov):
	'''Get the normalized PDF from the Asimov dataset (nominalMC + data-driven prediction) for given region.
	Args:
	    region (str): name of region to create PDF for
	    asimov (TH2): asimov dataset
	Returns:
	    PDF (TH2): The normalized total background PDF for the given region
	'''
	PDF = asimov.Clone('PDF_{}'.format(region))
	# The PDF must be scaled by the integral of teh *Asimov* dataset, but this is *NOT* the number of
	# events that should be generated, since that is not representative of the number of events in the
	# actual dataset, which is what we are trying to model
	nEventsAsimov = PDF.Integral(1, PDF.GetNbinsX()+1, 1, PDF.GetNbinsY()+1)
	PDF.Scale(1./nEventsAsimov)
	PDF.SetDirectory(0)
	return PDF

    def GetCDF(self, region, PDF):
	nx = PDF.GetNbinsX()
	ny = PDF.GetNbinsY()
	CDF = ROOT.TH1F('CDF_{}'.format(region),"",nx*ny,0,nx*ny)
	CDF.SetDirectory(0)
	cumulativeBin = 0
	for i in range(1, nx+1):
	    for j in range(1, ny+1):
		cumulativeBin += 1
		PDFval = PDF.GetBinContent(i,j) + CDF.GetBinContent(cumulativeBin-1)
		CDF.SetBinContent(cumulativeBin, PDFval)
	# Check that the CDF is bounded within [0, 1.)
	if (CDF.GetMaximum() > 1.0):
	    raise ValueError('CDF maximum ({}) > 1.0 - something has gone wrong'.format(CDF.GetMaximum()))
	return CDF

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

    def GeneratePseudoData(self, region, PDF, CDF, nEvents, name):
	toy = self.template_hist.Clone(name)
	toy.SetDirectory(0)
	print('Generating {} events for toy data hist {}'.format(int(nEvents),name))
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
	print('\tMultiplying ({}) x ({}) = ({})'.format(h1.GetName(),h2.GetName(),name))
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
