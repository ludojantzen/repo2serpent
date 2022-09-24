###############################################################################

# Serpent 2 Makefile

###############################################################################

# GNU Compiler:

CC  	 = gcc
CFLAGS   = -Wall -ansi -ffast-math -O3 -Wunused
LDFLAGS  = -lm

# Define installation directory for make install:

#PREFIX = ...

# Parallel calculation using Open MP:

CFLAGS  += -DOPEN_MP
CFLAGS  += -fopenmp
LDFLAGS += -fopenmp

# This is needed in newer gcc versions to supress some unnecessary warnings

CFLAGS += -Wno-unused-but-set-variable

# Remove this if compilation with mpicc produces unnecessary warnings

CFLAGS += -pedantic

# Flags for shared library compilation

CFLAGS += -fpic

# Additional flags for development

#CFLAGS += -Wshadow
#CFLAGS += -Wconversion -Wno-sign-conversion
#CFLAGS += -Wconversion -Wno-sign-conversion -Wduplicated-cond -Wlogical-op -Wnull-dereference -Wjump-misses-init -Wdouble-promotion -Wformat=2
#CFLAGS += -Wextra

# This is needed by new version of XCode in OS/X

#CFLAGS  += -D_FORTIFY_SOURCE=0

###############################################################################

# Intel Compiler:

#CC  	 = icc

# Optimization for Intel CPU's:

#CFLAGS   = -Wall -ansi -pedantic -xHost -ipo -DINTELCC

# Alternative:

#CFLAGS   = -Wall -ansi -pedantic -O3 -DINTELCC

#LDFLAGS  = -lm

# Parallel calculation using Open MP:

#CFLAGS  += -DOPEN_MP
#CFLAGS  += -openmp
#LDFLAGS += -openmp

###############################################################################

# Parallel calculation using MPI:

# NOTE: The use of hybrid MPI/OpenMP mode requires thread-safe MPI
#       implementation. Some MPI implementations, such as some versions (?) of
#       Open MPI are not thread safe, which will cause problems in memory
#       management routines (calloc, realloc and free). These problems may
#       result in failure in memory allocation or unexpected behaviour due to
#       corrupted registers (?).

#CC  	 = mpicc
#CFLAGS  += -DMPI

###############################################################################

# GD graphics library:

LDFLAGS += -lgd

# Compile with the following option if GD library is not available:

#CFLAGS += -DNO_GFX_MODE

###############################################################################

# Debugging:

#CFLAGS  += -DDEBUG

#CFLAGS  += -g
#LDFLAGS += -g

# Profiler:

#CFLAGS  += -pg
#LDFLAGS += -pg

###############################################################################

OBJS	=	addbranching.o \
		addbuf.o \
		addbuf1d.o \
		addchains.o \
		addfet.o \
		addddres.o \
		additem.o \
		addmesh.o \
		addmeshidx.o \
		addnuclide.o \
		addprivatedata.o \
		addprivateres.o \
		addpts.o \
		addsabdata.o \
		addsearchmesh.o \
		addsortitem.o \
		addstablenuclides.o \
		addstat.o \
		addstlpoint.o \
		addvaluepair.o \
		adjustenergygrid.o \
		adjustrmx.o \
		adjustsabdata.o \
		allocfetcache.o \
		allocinterfacestat.o \
		allocmacroxs.o \
		allocmicroxs.o \
		allocparticlestack.o \
		allocprecdet.o \
		allocprivatedata.o \
		allocstathistory.o \
		allocvaluepair.o \
		alpha.o \
		alphaxs.o \
		applygcsymmetries.o \
		atof.o \
		atoi.o \
		atomicrelaxation.o \
		averagetransmuxs.o \
		avgrmx.o \
		azirot.o \
		b1flux.o \
		b1fluxcorr.o \
		bankstostore.o \
		branchfrac.o \
		broadcastifcdata.o \
		broadcrosssection.o \
		boundaryconditions.o \
		bsfunn.o \
		bufmean.o \
		bufn.o \
		bufval.o \
		bufwgt.o \
		burnmatcompositions.o \
		burnmaterials.o \
		burnmatrixsize.o \
		burnupcycle.o \
		cachexs.o \
		calcconvcriteria.o \
		calcmicrogroupxs.o \
		calculateactivities.o \
		calculatebytes.o \
		calculatedtmajorants.o \
		calculateentropies.o \
		calculatelegendremoments.o \
		calculatemasses.o \
		calculatemgxs.o \
		calculaterelalpha.o \
		calculaterelpopsize.o \
		calculatetetboundingbox.o \
		calculatetetcenter.o \
		calculatetransmuxs.o \
		calculateuresmajorants.o \
		ccsmatrixcolspace.o \
		ccsmatrixcopy.o \
		ccsmatrixfree.o \
		ccsmatrixisort.o \
		ccsmatrixnew.o \
		ccsmatrixprint.o \
		cellcount.o \
		cellvolumes.o \
		checkcoefcalc.o \
		checkdddomain.o \
		checkddrecv.o \
		checkduplicates.o \
		checkfinishddasynch.o \
		checkfinishddsynch.o \
		checknuclidedata.o \
		checkpointer.o \
		checkpolyhedmesh.o \
		checkvalue.o \
		checkrealistsum.o \
		checkunused.o \
		cleanupddcomms.o \
		cleanupddrecv.o \
		cleanupddsend.o \
		clearbuf.o \
		clearinterfacestat.o \
		clearmicrogroupxs.o \
		clearprivatedata.o \
		clearprivateres.o \
		clearreltransmuxs.o \
		clearstat.o \
		cleartransmuxs.o \
		clearvaluepairs.o \
		closelist.o \
		coldet.o \
		collectbuf.o \
		collectburndata.o \
		collectdet.o \
		collectfet.o \
		collectdyndata.o \
		collectparalleldata.o \
		collectprecdet.o \
		collectresults.o \
		collectsensresults.o \
		collectuncresults.o \
		collectvrmeshdata.o \
		collision.o \
		compresssenslabel.o \
		coefcycle.o \
		coefoutput.o \
		combineactinides.o \
		combinefissionyields.o \
		commandlineplotter.o \
		comparestr.o \
		complex.o \
		complexrea.o \
		comptonscattering.o \
		contribsplit.o \
		coordexpans.o \
		coordtrans.o \
		countdynsrc.o \
		countsensparams.o \
		covmatrixfromblock.o \
		covmatrixfromsingle.o \
		createfinixifc.o \
		creategeometry.o \
		createmesh.o \
		createumshcells.o \
		createuniverse.o \
		cspline.o \
		cyldis.o \
		dataifcadens.o \
		dataifcxs.o \
		decaymeshprecdet.o \
		decaypointprecdet.o \
		decomposeelements.o \
		defaultbradata.o \
		deinitsocket.o \
		densityfactor.o \
		depletionpolyfit.o \
		detbin.o \
		detectoroutput.o \
		deterministicleakage.o \
		detidx.o \
		detresponse.o \
		dfpos.o \
		dfsol.o \
		dfsolver.o \
		die.o \
		diffcoefed.o \
		disperse.o \
		disperse2.o \
		distributeddmajorants.o \
		distributeddmaxadens.o \
		distributeddlimbo.o \
		distributematerialdata.o \
		divideburnmat.o \
		dividemeshcell.o \
		dividepolyhedcell.o \
		dividepolyhedface.o \
		dividezone.o \
		divsubdistance.o \
		dopmicroxs.o \
		dopplerbroad.o \
		dtmajorant.o \
		duplicateitem.o \
		duplicateparticle.o \
		eblockfrombank.o \
		eblocktobank.o \
		elasticscattering.o \
		elbrsxsdata.o \
		element.o \
		elspdata.o \
		endfcolf.o \
		endfcoli.o \
		endfinterp.o \
		endfmat.o \
		endfmf.o \
		endfmt.o \
		endfnewline.o \
		eneimportance.o \
		error.o \
		estimateruntime.o \
		eventfrombank.o \
		eventstosensitivity.o \
		eventtobank.o \
		expandfe.o \
		expandprivatearrays.o \
		expodecfit.o \
		extractsenslabel.o \
		feifc.o \
		fetidx.o \
		fetfinalize.o \
		fillstlmesh.o \
		finalizeccprecdet.o \
		finalizempi.o \
		finalizeumsh.o \
		findregmeshifcindex.o \
		findsenseindex.o \
		findsensindices.o \
		findsensmuindex.o \
		findsenszaiindex.o \
		findtetcell.o \
		findinterfaceregions.o \
		findlatticeregion.o \
		findmaterialpointers.o \
		findnestregion.o \
		findnuclidedata.o \
		findpbregion.o \
		findrowindexes.o \
		findstlsolid.o \
		finduniversecell.o \
		findxslimits.o \
		finixptrfromfpe.o \
		firstitem.o \
		fisse.o \
		fissecomp.o \
		fissmtxindex.o \
		fissmtxoutput.o \
		fission.o \
		fixhexmesh.o \
		fixpolyhedmesh.o \
		flushbank.o \
		flushprecsource.o \
		flushddparticles.o \
		formtransmupaths.o \
		freeddcomms.o \
		freefinix.o \
		freemem.o \
		frombank.o \
		fromcommonque.o \
		fromsrc.o \
		fromstack.o \
		fromstore.o \
		fromtrkbank.o \
		fromque.o \
		gaussiansubst.o \
		geoimportance.o \
		geometryplotter.o \
		getbankedprecursors.o \
		getburnids.o \
		getifpancestordata.o \
		getrequiredcovmatrices.o \
		getsenseventtype.o \
		getimportantpts.o \
		getlatticeindexes.o \
		getlinenumber.o \
		getparams.o \
		getprivatedata.o \
		getprivateres.o \
		getsignalfromsupervisor.o \
		gettemp.o \
		gettext.o \
		gridfactor.o \
		gridsearch.o \
		gsconfigdata.o \
		hessfactorization.o \
		hexrotatecell.o \
		hexrotateface.o \
		homoflux.o \
		hismean.o \
		hisrelerr.o \
		hisval.o \
		icmidx.o \
		idxstr.o \
		ifcpoint.o \
		importancesolver.o \
		incell.o \
		inelasticscattering.o \
		initdata.o \
		initddcomms.o \
		initddrecv.o \
		inithistories.o \
		initialcritsrc.o \
		initmpi.o \
		initomp.o \
		initprecdet.o \
		initprecdetsource.o \
		initsignal.o \
		initsocket.o \
		insupercell.o \
		interpolatedata.o \
		interpolatefisse.o \
		interpolatemupdf.o \
		interpolatenubar.o \
		interpolatesab.o \
		intersectionlist.o \
		intetcell.o \
		inumshcell.o \
		invokebranch.o \
		involute.o \
		isotopefractions.o \
		isotozai.o \
		isotropicdirection.o \
		iteratecc.o \
		iterateexternal.o \
		iteratefinix.o \
		iteratekeff.o \
		iternucxs.o \
		kleinnishina.o \
		lagrangeinterp.o \
		lastitem.o \
		latcellorig.o \
		leak.o \
		leakdet.o \
		levelscattering.o \
		lifolistsize.o \
		linkgcumaterials.o \
		linkinterfacematerials.o \
		linkreactions.o \
		linksabdata.o \
		listptr.o \
		ludecomposition.o \
		main.o \
		majorantxs.o \
		makearray.o \
		makeburnmatrix.o \
		makeburnmatrixmsr.o \
		makedepletionzones.o \
		makeenergygrid.o \
		makepalette.o \
		makering.o \
		materialburnup.o \
		materialtotals.o \
		materialvolumes.o \
		matlaboutput.o \
		matproduct.o \
		matptr.o \
		matpos.o \
		matrixexponential.o \
		maxsrcimp.o \
		maxsurfdimensions.o \
		maxwellenergy.o \
		mean.o \
		mem.o \
		memcount.o \
		meshcellbounds.o \
		meshcellconnectdiags.o \
		meshcellfromcgns.o \
		meshcellgetface.o \
		meshcellgetfacepos.o \
		meshcellindirection.o \
		meshcellrotatelists.o \
		meshcelltype.o \
		meshcellvol.o \
		meshindex.o \
		meshplotter.o \
		meshptr.o \
		meshtot.o \
		meshval.o \
		macrourescorr.o \
		macroxs.o \
		mgxs.o \
		microcalc.o \
		microdepoutput.o \
		micromajorantxs.o \
		microxs.o \
		minxs.o \
		moraoutput.o \
		movedt.o \
		moveifc.o \
		moveitemfirst.o \
		moveitemlast.o \
		moveitemright.o \
		movest.o \
		movestore.o \
		mpitransfer.o \
		msrrealist.o \
		myparallelmat.o \
		nearestboundary.o \
		nearestpbsurf.o \
		nearestmeshboundary.o \
		neareststlsurf.o \
		nearestumshsurf.o \
		nestvolumes.o \
		newitem.o \
		newlifoitem.o \
		newrealist.o \
		newstat.o \
		nextitem.o \
		nextreaction.o \
		nextword.o \
		nfkerma.o \
		nocorrector.o \
		normcoef.o \
		normalizecompositions.o \
		normalizecritsrc.o \
		normalizedynsrc.o \
		normalizeimp.o \
		normalizeprecdet.o \
		note.o \
		nubar.o \
		numericgauss.o \
		numericstr.o \
		nxn.o \
		ompresetcomp.o \
		ompsetcomp.o \
		omptestcomp.o \
		opendatafile.o \
		otfburntta.o \
		otfburnxs.o \
		otfsabxs.o \
		otfsabscattering.o \
		overrideids.o \
		particlesfromstore.o \
		pairproduction.o \
		parlett.o \
		parsecommandline.o \
		pause.o \
		phidis.o \
		photoelectric.o \
		photonmacroxs.o \
		photonmicroxs.o \
		photonprod.o \
		plotimage.o \
		plottracks.o \
		poisoneq.o \
		poisonxs.o \
		polarangle.o \
		polynomiallegendre.o \
		polynomialzernike.o \
		polypinf.o \
		polysameface.o \
		posannihilation.o \
		postddrecv.o \
		postddsend.o \
		potcorr.o \
		preallocmem.o \
		precdet.o \
		precursorpopcontrol.o \
		precursorstostore.o \
		preparecciter.o \
		preparecellsearchmesh.o \
		preparetransportcycle.o \
		presort.o \
		pretrans.o \
		previtem.o \
		printcellstats.o \
		printcoevals.o \
		printcompositions.o \
		printcoredistr.o \
		printcycleoutput.o \
		printdepmatrix.o \
		printdepoutput.o \
		printdepvals.o \
		printfinix.o\
		printgammaspectra.o \
		printgeometrydata.o \
		printhistoryoutput.o \
		printinterfaceoutput.o \
		printmaterialdata.o \
		printmeshcell.o \
		printmixtures.o \
		printmvar.o \
		printnuclidedata.o \
		printpbdata.o \
		printprecdet.o \
		printprogress.o \
		printsrcimp.o \
		printreactionlists.o \
		printtmsdiagnostics.o \
		printtitle.o \
		printvalues.o \
		probevrmesh.o \
		processbc.o \
		processbuinterface.o \
		processburnmat.o \
		processburnupegroups.o \
		processcells.o \
		processcellmesh.o \
		processcovariancedata.o \
		processcpd.o \
		processcomplementcells.o \
		processcompton.o \
		processdatainterfaces.o \
		processdecaydata.o \
		processdecaysrc.o \
		processdephis.o \
		processdetectors.o \
		processdivisors.o \
		processdix.o \
		processelectrons.o \
		processentropy.o \
		processedistributions.o \
		processevents.o \
		processfinix.o \
		processfissmtx.o \
		processfissionyields.o \
		processgc.o \
		processicm.o \
		processifcfb.o \
		processifcfetmaterial.o \
		processifcfunc.o \
		processifcptavg.o \
		processifcregmesh.o \
		processifctetmesh.o \
		processinelaxs.o \
		processinterface.o \
		processinventory.o \
		processiternucs.o \
		processlattices.o \
		processmaterials.o \
		processmeshplots.o \
		processmixture.o \
		processmsr.o \
		processmudistributions.o \
		processnests.o \
		processnubardata.o \
		processnuclides.o \
		processotfburn.o \
		processpairproduction.o \
		processphotoelectric.o \
		processpbgeometry.o \
		processphotonatt.o \
		processphotonecut.o \
		processphotonprod.o \
		processphotonrea.o \
		processpoisons.o \
		processprecdet.o \
		processrayleigh.o \
		processreactionlists.o \
		processrelaxation.o \
		processreprocessors.o \
		processrmx.o \
		processsenseblocks.o \
		processsensitivities.o \
		processsensstats.o \
		processsingleinterface.o \
		processsources.o \
		processstats.o \
		processstlgeometry.o \
		processsurfaces.o \
		processsymmetries.o \
		processtmpdata.o \
		processtimebins.o \
		processtransformations.o \
		processtranspcorr.o \
		processttb.o \
		processumshgeometry.o \
		processuresdata.o \
		processuseregrids.o \
		processvr.o \
		processxsdata.o \
		pulsedet.o \
		putprivatedata.o \
		putcompositions.o \
		putmeshidx.o \
		putotfburnconc.o \
		putplotcolors.o \
		putpoisonconc.o \
		puttext.o \
		pythonplotter.o \
		qrfactorization.o \
		radgammasrc.o \
		rand64.o \
		randf.o \
		randnorm.o \
		rayleighscattering.o \
		reactioncount.o \
		reactioncutoff.o \
		reactionmt.o \
		reactiontargetzai.o \
		readacefile.o \
		readbrafile.o \
		readcoverxfile.o \
		readdatainterfaces.o \
		readdecayfile.o \
		readdirectoryfile.o \
		readfinixifc.o \
		readfissionyields.o \
		readifcbins.o \
		readifcfb.o \
		readifcfblims.o \
		readifcfe.o \
		readifcfunc.o \
		readifcofmesh.o \
		readifcptavg.o \
		readifcregmesh.o \
		readifctetmesh.o \
		readinfix.o \
		readinput.o \
		readinterface.o \
		readmesh.o \
		readmeshptr.o \
		readpendfdata.o \
		readofpatches.o \
		readofdata.o \
		readofdensities.o \
		readoffaces.o \
		readofheader.o \
		readofmapping.o \
		readofmaterials.o \
		readofmesh.o \
		readofownnbr.o \
		readofpoints.o \
		readoftemperatures.o \
		readpbgeometry.o \
		readplasmasrc.o \
		readphotondata.o \
		readrestartfile.o \
		readsourcefile.o \
		readstlgeometry.o \
		readstlmesh.o \
		readtextfile.o \
		readumshgeometry.o \
		readwwdmesh.o \
		reallocmem.o \
		reallocddsend.o \
		reamulti.o \
		receiveifcinputdata.o \
		receiveddparticlecounters.o \
		receiveddparticles.o \
		recoildet.o \
		recyclermxdata.o \
		redistributeques.o \
		redistributestacks.o \
		reducebuffer.o \
		reduceprivateres.o \
		refreshinventory.o \
		reinitrng.o \
		relaxccresults.o \
		relaxinterfacepower.o \
		relaxtransmuxs.o \
		relerr.o \
		removeflaggeditems.o \
		removeitem.o \
		removestorefiles.o \
		removevoidcells.o \
		reopenlist.o \
		replaceitem.o \
		replacephotondata.o \
		reprocess.o \
		resetoption.o \
		resetotfburn.o \
		resetpoisonconc.o \
		resettimer.o \
		resetddcomms.o \
		resetddsend.o \
		resizedynsrc.o \
		resizefissionsrc.o \
		responsefunction.o \
		riacycle.o \
		rroutput.o \
		runfinix.o \
		sabscattering.o \
		sampledelnu.o \
		sampleendflaw.o \
		sampleifcdata.o \
		samplemeshdelnu.o \
		samplemu.o \
		samplenu.o \
		sampleplasmasrc.o \
		samplepointdelnu.o \
		sampleprecursorgroup.o \
		sampleptable.o \
		samplereaction.o \
		samplesrcpoint.o \
		sampletabular.o \
		scalarprod.o \
		schurfactorization.o \
		score.o \
		scorealb.o \
		scoreactidet.o \
		scoredf.o \
		scorecapture.o \
		scorecpd.o \
		scorecmm.o \
		scoredirectsensitivity.o \
		scorefet.o \
		scorefission.o \
		scoregc.o \
		scoreicmcol.o \
		scoreicmtrk.o \
		scoreinterfaceflux.o \
		scoreinterfacepower.o \
		scoreiternucs.o \
		scoremesh.o \
		scoremicrodep.o \
		scorermxcurr.o \
		scorermxmfp.o \
		scorermxresp.o \
		scorermxsrc.o \
		scoreotfburn.o \
		scorepb.o \
		scoreperturbations.o \
		scorephotonheat.o \
		scorepinpower.o \
		scorepoison.o \
		scorescattering.o \
		scorescaresp.o \
		scoresenslabel.o \
		scoresurf.o \
		scoretimeconstants.o \
		scoretimesource.o \
		scoretransmuxs.o \
		scoreufs.o \
		searcharray.o \
		seekarray.o \
		seeklist.o \
		seekliststr.o \
		sendddparticle.o \
		sendddparticlecounter.o \
		sendifcinputtemplates.o \
		sendifcmesh.o \
		sendifcoutputdata.o \
		sendifcoutputtemplates.o \
		sensbufval.o \
		senscollision.o \
		sensitivityoutput.o \
		separatebanuc.o \
		setcipopulation.o \
		setcoefcalc.o \
		setddidsimple.o \
		setdepstepsize.o \
		setdetflags.o \
		setdirectpointers.o \
		setfisse.o \
		sethisvcalc.o \
		setifctmslimits.o \
		setnormalization.o \
		setoption.o \
		setoptimization.o \
		setpathlevels.o \
		setprecursorgroups.o \
		setstlmeshpointers.o \
		shareinputdata.o \
		shuntingyard.o \
		signalexternal.o \
		signalhandler.o \
		socketreceive.o \
		socketreceivedouble.o \
		socketreceivelong.o \
		socketreceivestring.o \
		socketsend.o \
		socketsenddouble.o \
		socketsendlong.o \
		socketsendstring.o \
		solveotfburn.o \
		solvermx.o \
		solvermxcrit.o \
		sortall.o \
		sortarray.o \
		sortlist.o \
		speed.o \
		splitlist.o \
		srcdet.o \
		srcroulette.o \
		starttimer.o \
		statbin.o \
		statsum.o \
		stattests.o \
		stdcomp.o \
		stddev.o \
		stlfacetdistance.o \
		stlmatfinder.o \
		stlraytest.o \
		stopatboundary.o \
		stopatwwbound.o \
		stopci.o \
		stopcciter.o \
		stopsurfdetflg.o \
		stoptimer.o \
		storecomposition.o \
		storehistorypoint.o \
		storermxevent.o \
		storesensevent.o \
		storetransmuxs.o \
		storevaluepair.o \
		sumdivcompositions.o \
		sumprivatedata.o \
		sumprivateres.o \
		sumtotxs.o \
		superdet.o \
		supervisedcycle.o \
		surfacedistance.o \
		surfacenormal.o \
		surfacesrc.o \
		surfacevol.o \
		swapitems.o \
		swapuniverses.o \
		symboliclu.o \
		symmetryboundary.o \
		systemstat.o \
		targetvelocity.o \
		tmpmajorants.o \
		testasciifile.o \
		testdosfile.o \
		testfacetoverlap.o \
		testhisvbreak.o \
		testparam.o \
		testsurface.o \
		teststlgeometry.o \
		teststlsolids.o \
		testunisym.o \
		testvaluepair.o \
		testxs.o \
		tetputboundingbox.o \
		tetravol.o \
		thingrid.o \
		timecutoff.o \
		timeintervalstr.o \
		timercpuval.o \
		timerval.o \
		timestamp.o \
		timestr.o \
		tobank.o \
		tocommonque.o \
		tolimbo.o \
		torusdis.o \
		toque.o \
		tostack.o \
		tostore.o \
		totxs.o \
		trackfile.o \
		tracking.o \
		trackingerror.o \
		trackmode.o \
		transportcycle.o \
		transportcorrection.o \
		trapz.o \
		trapzreal.o \
		truncate.o \
		tta.o \
		ttb.o \
		ttachain.o \
		ttaloop.o \
		ufsfactor.o \
		unionizegrid.o \
		unisym.o \
		universeboundaries.o \
		updatecistop.o \
		updatefinixifc.o \
		updatefinixpower.o \
		updateifcdensmax.o \
		updateifctempminmax.o \
		updatemicrodens.o \
		updatermxwgt.o \
		uresdilumicroxs.o \
		uresfactor.o \
		userifc.o \
		useriter.o \
		usersrc.o \
		usersurf.o \
		usertransmucut.o \
		valuepairidx.o \
		valuepairval.o \
		vectornorm.o \
		vectorprod.o \
		virtgcucolflags.o \
		volumesmc.o \
		vrcycle.o \
		walkeralias.o \
		warn.o \
		weightwindow.o \
		whereami.o \
		workarray.o \
		writecimomfluxes.o \
		writedepfile.o \
		writedynsrc.o \
		writefinixinputfile.o \
		writetetmeshtogeo.o \
		writefinixifc.o \
		writeicmdata.o \
		writesensresult.o \
		writesourcefile.o \
		writestlmesh.o \
		writeumshtostl.o \
		writewwmesh.o \
		wwdis.o \
		wwimportance.o \
		wwinsrc.o \
		xsplotter.o \
		zaitoiso.o \
		zdis.o \
		zonecount.o

FOBJS	=	aux_functions.o \
		clmech.o  \
		clmechaux.o \
		clmechprop.o \
		clthprop.o \
		coolant.o \
		finixdata.o \
		finix_initialization.o \
		finix_output.o \
		fumech.o \
		fumechprop.o \
		futhprop.o \
		gap.o \
		heateq1d.o \
		initial.o \
		transient.o

OBSOLETE = 	capture.c \
		nextlistrea.c \
		scoreadjoint.c \
		statn.c \
		rewindrealist.c \
		contribdet.c \
		uresdilutioncutoff.c \
		currentdet.c \
		rendezvous.c \
		scoredet.c \
		tracktimecut.c \
		checkspecialinventoryentry.c \
		privatedataptr.c \
		preallocxsmem.c \
		tosrc.c \
		transferparticle.c \
		initpartstacks.c \
		readtransmudata.c \
		generatetransmuchain.c \
		processontheflydata.c \
		processdbrcdata.c \
		findtargetnuclides.c \
		p1solver.c \
		scoreadf.c \
		userdf.c \
		calculatemajorant.c \
		flushsrc.c \
		readibrfile.c \
		albedoidx.c \
		albedooutput.c \
		processfum.c \
		scoreicm.c \
		addtrackpt.c \
		storetimepts.c \
		storecompositions.c \
		testplane.c \
		sampletarget.c \
		cgnscellcenterpoint.c \
		dividecgnscell.c \
		dividecgnsface.c \
		findcgnscell.c \
		fixcgnsmesh.c \
		writecgnstogeo.c \
		collectfinix.c \
		hdmc.c \
		clearevents.c \
		thompsonscattering.c \
		sabdatatester.c \
		polyhedcellcenterpoint.c \
		updateinterface.c \
		profile.c \
		ucbburnupcycle.c \
		collecttfb.c \
		identifytfbmat.c \
		iteratetfb.c \
		printtfboutput.c \
		processtfb.c \
		reinittfb.c \
		scoretemp.c \
		updatetfbgeom.c \
		updatetfbhc.c \
		hexfromcgns.c \
		dividehexcell.c \
		hexgetface.c \
		hexgetfacepos.c \
		hexconnectdiags.c \
		hexindirection.c \
		hexrotatelists.c \
                database.c \
		db_functions.c \
		defaults.c \
		finixtextfile.c \
		nextline.c \
		readfile.c \
		caxoutput.c \
		storesimdata.c \
		scorewwdcurr.c \
		processbradata.c \
		distributefinix.c \
		vrtester.c \
		aresoutput.c \
		hexnewtetface.c \
		readmaterialcomp.c \
		retrievecomposition.c \
		b1solver.c \
		distributelimbo.c \
		setdecompid.c \
		sddcomms.c \
		simplifyumshsurfaces.c \
		readofbatches.c \
		collectmicrodepdata.o \

###############################################################################

# Compilation with FINIX:

#OBJS += $(FOBJS)
#CFLAGS  += -DFINIX

###############################################################################

EXE     = sss2
RM      = /bin/rm -f
ZIP	= zip
TAR	= tar -zcf

DATE 	= -DMC_COMPILE_DATE="\"`date `\""

# Version number added:

BKFILE 	= "`date +'backup/serpent_2.1.32_%y%m%d%H%M.zip'`"
EXFILE  = "`date +'backup/serpent_2.1.32_%y%m%d%H%M'`"

# Get number of bits:

shell = $SHELL
lbits := $(shell getconf LONG_BIT)

# Use this to override the 64 bit check

#lbits = 64

###############################################################################

all:    check $(EXE)
	@echo "Serpent 2 Compiled OK."

$(EXE): $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS) -o $(EXE)

lib:	$(OBJS)
	$(CC) -shared -o lib$(EXE).so $(OBJS) $(LDFLAGS)

$(DBG): $(OBJS)
	$(CC) -DDEBUG $(OBJS) $(LDFLAGS) -o $(EXE)

bk:
	$(ZIP) $(BKFILE) *.c *.h Makefile versions.txt
	cp $(EXE) $(EXFILE)
	cp $(BKFILE) ./serpent2.zip
	chmod a-w $(BKFILE)
	chmod a-w $(EXFILE)

tgz:
	$(TAR) serpent2.tar.gz *.c *.h Makefile versions.txt $(EXE)

clean:
	$(RM) $(OBJS) $(FOBJS) $(OBSOLETE)

# Check that system is 64-bit. NOTE: If this check fails in a 64-bit system,
# remove the exit command and report to: Jaakko.Leppanen@vtt.fi

check:

ifneq ($(lbits), 64)
	@echo "Serpent 2 must be compiled in a 64-bit system"
	@exit 2
endif

# Installation if installation directory is set

install: $(EXE)

ifeq ($(PREFIX),)
	@echo "No installation directory given. Modify Makefile"
	@exit 2
else
	@echo "Installing Serpent 2 to $(PREFIX)"
	install -m 750 -d $(PREFIX)/
	install -m 750 -d $(PREFIX)/bin/
	install -m 750 $(EXE) $(DESTDIR)$(PREFIX)/bin/
endif

###############################################################################

addbranching.o: addbranching.c header.h locations.h
	$(CC) $(CFLAGS) -c addbranching.c

addbuf.o: addbuf.c header.h locations.h
	$(CC) $(CFLAGS) -c addbuf.c

addbuf1d.o: addbuf1d.c header.h locations.h
	$(CC) $(CFLAGS) -c addbuf1d.c

addchains.o: addchains.c header.h locations.h
	$(CC) $(CFLAGS) -c addchains.c

addddres.o: addddres.c header.h locations.h
	$(CC) $(CFLAGS) -c addddres.c

additem.o: additem.c header.h locations.h
	$(CC) $(CFLAGS) -c additem.c

addfet.o: addfet.c header.h locations.h
	$(CC) $(CFLAGS) -c addfet.c

addmesh.o: addmesh.c header.h locations.h
	$(CC) $(CFLAGS) -c addmesh.c

addmeshidx.o: addmeshidx.c header.h locations.h
	$(CC) $(CFLAGS) -c addmeshidx.c

addnuclide.o: addnuclide.c header.h locations.h
	$(CC) $(CFLAGS) -c addnuclide.c

addprivatedata.o: addprivatedata.c header.h locations.h
	$(CC) $(CFLAGS) -c addprivatedata.c

addprivateres.o: addprivateres.c header.h locations.h
	$(CC) $(CFLAGS) -c addprivateres.c

addpts.o: addpts.c header.h
	$(CC) $(CFLAGS) -c addpts.c

addsabdata.o: addsabdata.c header.h locations.h
	$(CC) $(CFLAGS) -c addsabdata.c

addsearchmesh.o: addsearchmesh.c header.h locations.h
	$(CC) $(CFLAGS) -c addsearchmesh.c

addsortitem.o: addsortitem.c header.h locations.h
	$(CC) $(CFLAGS) -c addsortitem.c

addstablenuclides.o: addstablenuclides.c header.h locations.h nuclide_masses.h
	$(CC) $(CFLAGS) -c addstablenuclides.c

addstat.o: addstat.c header.h locations.h
	$(CC) $(CFLAGS) -c addstat.c

addstlpoint.o: addstlpoint.c header.h locations.h
	$(CC) $(CFLAGS) -c addstlpoint.c

addvaluepair.o: addvaluepair.c header.h locations.h
	$(CC) $(CFLAGS) -c addvaluepair.c

adjustenergygrid.o: adjustenergygrid.c header.h locations.h
	$(CC) $(CFLAGS) -c adjustenergygrid.c

adjustrmx.o: adjustrmx.c header.h locations.h
	$(CC) $(CFLAGS) -c adjustrmx.c

adjustsabdata.o: adjustsabdata.c header.h locations.h
	$(CC) $(CFLAGS) -c adjustsabdata.c

allocfetcache.o: allocfetcache.c header.h locations.h
	$(CC) $(CFLAGS) -c allocfetcache.c

allocinterfacestat.o: allocinterfacestat.c header.h locations.h
	$(CC) $(CFLAGS) -c allocinterfacestat.c

allocmacroxs.o: allocmacroxs.c header.h locations.h
	$(CC) $(CFLAGS) -c allocmacroxs.c

allocmicroxs.o: allocmicroxs.c header.h locations.h
	$(CC) $(CFLAGS) -c allocmicroxs.c

allocparticlestack.o: allocparticlestack.c header.h locations.h
	$(CC) $(CFLAGS) -c allocparticlestack.c

allocprecdet.o: allocprecdet.c header.h locations.h
	$(CC) $(CFLAGS) -c allocprecdet.c

allocprivatedata.o: allocprivatedata.c header.h locations.h
	$(CC) $(CFLAGS) -c allocprivatedata.c

allocstathistory.o: allocstathistory.c header.h locations.h
	$(CC) $(CFLAGS) -c allocstathistory.c

allocvaluepair.o: allocvaluepair.c header.h locations.h
	$(CC) $(CFLAGS) -c allocvaluepair.c

alpha.o: alpha.c header.h locations.h
	$(CC) $(CFLAGS) -c alpha.c

alphaxs.o: alphaxs.c header.h locations.h
	$(CC) $(CFLAGS) -c alphaxs.c

applygcsymmetries.o: applygcsymmetries.c header.h locations.h
	$(CC) $(CFLAGS) -c applygcsymmetries.c

atof.o: atof.c header.h
	$(CC) $(CFLAGS) -c atof.c

atoi.o: atoi.c header.h
	$(CC) $(CFLAGS) -c atoi.c

atomicrelaxation.o: atomicrelaxation.c header.h locations.h
	$(CC) $(CFLAGS) -c atomicrelaxation.c

averagetransmuxs.o: averagetransmuxs.c header.h locations.h
	$(CC) $(CFLAGS) -c averagetransmuxs.c

avgrmx.o: avgrmx.c header.h locations.h
	$(CC) $(CFLAGS) -c avgrmx.c

azirot.o: azirot.c header.h
	$(CC) $(CFLAGS) -c azirot.c

b1flux.o: b1flux.c header.h locations.h
	$(CC) $(CFLAGS) -c b1flux.c

b1fluxcorr.o: b1fluxcorr.c header.h locations.h
	$(CC) $(CFLAGS) -c b1fluxcorr.c

bankstostore.o: bankstostore.c header.h locations.h
	$(CC) $(CFLAGS) -c bankstostore.c

branchfrac.o: branchfrac.c header.h locations.h
	$(CC) $(CFLAGS) -c branchfrac.c

broadcastifcdata.o: broadcastifcdata.c header.h locations.h
	$(CC) $(CFLAGS) -c broadcastifcdata.c

broadcrosssection.o: broadcrosssection.c header.h locations.h
	$(CC) $(CFLAGS) -c broadcrosssection.c

boundaryconditions.o: boundaryconditions.c header.h locations.h
	$(CC) $(CFLAGS) -c boundaryconditions.c

bsfunn.o: bsfunn.c header.h locations.h
	$(CC) $(CFLAGS) -c bsfunn.c

bufmean.o: bufmean.c header.h locations.h
	$(CC) $(CFLAGS) -c bufmean.c

bufn.o: bufn.c header.h locations.h
	$(CC) $(CFLAGS) -c bufn.c

bufval.o: bufval.c header.h locations.h
	$(CC) $(CFLAGS) -c bufval.c

bufwgt.o: bufwgt.c header.h locations.h
	$(CC) $(CFLAGS) -c bufwgt.c

burnmatcompositions.o: burnmatcompositions.c header.h locations.h
	$(CC) $(CFLAGS) -c burnmatcompositions.c

burnmaterials.o: burnmaterials.c header.h locations.h
	$(CC) $(CFLAGS) -c burnmaterials.c

burnmatrixsize.o: burnmatrixsize.c header.h locations.h
	$(CC) $(CFLAGS) -c burnmatrixsize.c

burnupcycle.o: burnupcycle.c header.h locations.h
	$(CC) $(CFLAGS) -c burnupcycle.c

cachexs.o: cachexs.c header.h locations.h
	$(CC) $(CFLAGS) -c cachexs.c

calcconvcriteria.o: calcconvcriteria.c header.h locations.h
	$(CC) $(CFLAGS) -c calcconvcriteria.c

calcmicrogroupxs.o: calcmicrogroupxs.c header.h locations.h
	$(CC) $(CFLAGS) -c calcmicrogroupxs.c

calculateactivities.o: calculateactivities.c header.h locations.h
	$(CC) $(CFLAGS) -c calculateactivities.c

calculatebytes.o: calculatebytes.c header.h locations.h
	$(CC) $(CFLAGS) -c calculatebytes.c

calculatedtmajorants.o: calculatedtmajorants.c header.h locations.h
	$(CC) $(CFLAGS) -c calculatedtmajorants.c

calculateentropies.o: calculateentropies.c header.h locations.h
	$(CC) $(CFLAGS) -c calculateentropies.c

calculatelegendremoments.o: calculatelegendremoments.c header.h locations.h
	$(CC) $(CFLAGS) -c calculatelegendremoments.c

calculatemasses.o: calculatemasses.c header.h locations.h
	$(CC) $(CFLAGS) -c calculatemasses.c

calculatemgxs.o: calculatemgxs.c header.h locations.h
	$(CC) $(CFLAGS) -c calculatemgxs.c

calculaterelalpha.o: calculaterelalpha.c header.h locations.h
	$(CC) $(CFLAGS) -c calculaterelalpha.c

calculaterelpopsize.o: calculaterelpopsize.c header.h locations.h
	$(CC) $(CFLAGS) -c calculaterelpopsize.c

calculatetetboundingbox.o: calculatetetboundingbox.c header.h locations.h
	$(CC) $(CFLAGS) -c calculatetetboundingbox.c

calculatetetcenter.o: calculatetetcenter.c header.h locations.h
	$(CC) $(CFLAGS) -c calculatetetcenter.c

calculatetransmuxs.o: calculatetransmuxs.c header.h locations.h
	$(CC) $(CFLAGS) -c calculatetransmuxs.c

calculateuresmajorants.o: calculateuresmajorants.c header.h locations.h
	$(CC) $(CFLAGS) -c calculateuresmajorants.c

ccsmatrixcolspace.o: ccsmatrixcolspace.c header.h
	$(CC) $(CFLAGS) -c ccsmatrixcolspace.c

ccsmatrixcopy.o: ccsmatrixcopy.c header.h
	$(CC) $(CFLAGS) -c ccsmatrixcopy.c

ccsmatrixfree.o: ccsmatrixfree.c header.h
	$(CC) $(CFLAGS) -c ccsmatrixfree.c

ccsmatrixisort.o: ccsmatrixisort.c header.h
	$(CC) $(CFLAGS) -c ccsmatrixisort.c

ccsmatrixnew.o: ccsmatrixnew.c header.h
	$(CC) $(CFLAGS) -c ccsmatrixnew.c

ccsmatrixprint.o: ccsmatrixprint.c header.h
	$(CC) $(CFLAGS) -c ccsmatrixprint.c

cellcount.o: cellcount.c header.h locations.h
	$(CC) $(CFLAGS) -c cellcount.c

cellvolumes.o: cellvolumes.c header.h locations.h
	$(CC) $(CFLAGS) -c cellvolumes.c

checkcoefcalc.o: checkcoefcalc.c header.h locations.h
	$(CC) $(CFLAGS) -c checkcoefcalc.c

checkdddomain.o: checkdddomain.c header.h locations.h
	$(CC) $(CFLAGS) -c checkdddomain.c

checkddrecv.o: checkddrecv.c header.h locations.h
	$(CC) $(CFLAGS) -c checkddrecv.c

checkduplicates.o: checkduplicates.c header.h locations.h
	$(CC) $(CFLAGS) -c checkduplicates.c

checkfinishddasynch.o: checkfinishddasynch.c header.h locations.h
	$(CC) $(CFLAGS) -c checkfinishddasynch.c

checkfinishddsynch.o: checkfinishddsynch.c header.h locations.h
	$(CC) $(CFLAGS) -c checkfinishddsynch.c

checknuclidedata.o: checknuclidedata.c header.h locations.h
	$(CC) $(CFLAGS) -c checknuclidedata.c

checkpointer.o: checkpointer.c header.h locations.h
	$(CC) $(CFLAGS) -c checkpointer.c

checkpolyhedmesh.o: checkpolyhedmesh.c header.h locations.h
	$(CC) $(CFLAGS) -c checkpolyhedmesh.c

checkvalue.o: checkvalue.c header.h
	$(CC) $(CFLAGS) -c checkvalue.c

checkrealistsum.o: checkrealistsum.c header.h locations.h
	$(CC) $(CFLAGS) -c checkrealistsum.c

checkunused.o: checkunused.c header.h locations.h
	$(CC) $(CFLAGS) -c checkunused.c

cleanupddcomms.o: cleanupddcomms.c header.h locations.h
	$(CC) $(CFLAGS) -c cleanupddcomms.c

cleanupddrecv.o: cleanupddrecv.c header.h locations.h
	$(CC) $(CFLAGS) -c cleanupddrecv.c

cleanupddsend.o: cleanupddsend.c header.h locations.h
	$(CC) $(CFLAGS) -c cleanupddsend.c

clearbuf.o: clearbuf.c header.h locations.h
	$(CC) $(CFLAGS) -c clearbuf.c

clearinterfacestat.o: clearinterfacestat.c header.h locations.h
	$(CC) $(CFLAGS) -c clearinterfacestat.c

clearmicrogroupxs.o: clearmicrogroupxs.c header.h locations.h
	$(CC) $(CFLAGS) -c clearmicrogroupxs.c

clearprivatedata.o: clearprivatedata.c header.h locations.h
	$(CC) $(CFLAGS) -c clearprivatedata.c

clearprivateres.o: clearprivateres.c header.h locations.h
	$(CC) $(CFLAGS) -c clearprivateres.c

clearreltransmuxs.o: clearreltransmuxs.c header.h locations.h
	$(CC) $(CFLAGS) -c clearreltransmuxs.c

clearstat.o: clearstat.c header.h locations.h
	$(CC) $(CFLAGS) -c clearstat.c

cleartransmuxs.o: cleartransmuxs.c header.h locations.h
	$(CC) $(CFLAGS) -c cleartransmuxs.c

clearvaluepairs.o: clearvaluepairs.c header.h locations.h
	$(CC) $(CFLAGS) -c clearvaluepairs.c

closelist.o: closelist.c header.h locations.h
	$(CC) $(CFLAGS) -c closelist.c

coldet.o: coldet.c header.h locations.h
	$(CC) $(CFLAGS) -c coldet.c

collectbuf.o: collectbuf.c header.h locations.h
	$(CC) $(CFLAGS) -c collectbuf.c

collectburndata.o: collectburndata.c header.h locations.h
	$(CC) $(CFLAGS) -c collectburndata.c

collectdet.o: collectdet.c header.h locations.h
	$(CC) $(CFLAGS) -c collectdet.c

collectfet.o: collectfet.c header.h locations.h
	$(CC) $(CFLAGS) -c collectfet.c

collectdyndata.o: collectdyndata.c header.h locations.h
	$(CC) $(CFLAGS) -c collectdyndata.c

collectparalleldata.o: collectparalleldata.c header.h locations.h
	$(CC) $(CFLAGS) -c collectparalleldata.c

collectprecdet.o: collectprecdet.c header.h locations.h
	$(CC) $(CFLAGS) -c collectprecdet.c

collectresults.o: collectresults.c header.h locations.h
	$(CC) $(CFLAGS) -c collectresults.c

collectsensresults.o: collectsensresults.c header.h locations.h
	$(CC) $(CFLAGS) -c collectsensresults.c

collectuncresults.o: collectuncresults.c header.h locations.h
	$(CC) $(CFLAGS) -c collectuncresults.c

collectvrmeshdata.o: collectvrmeshdata.c header.h locations.h
	$(CC) $(CFLAGS) -c collectvrmeshdata.c

collision.o: collision.c header.h locations.h
	$(CC) $(CFLAGS) -c collision.c

compresssenslabel.o: compresssenslabel.c header.h locations.h
	$(CC) $(CFLAGS) -c compresssenslabel.c

comptonscattering.o: comptonscattering.c header.h locations.h
	$(CC) $(CFLAGS) -c comptonscattering.c

coefcycle.o: coefcycle.c header.h locations.h
	$(CC) $(CFLAGS) -c coefcycle.c

coefoutput.o: coefoutput.c header.h locations.h
	$(CC) $(CFLAGS) -c coefoutput.c

combineactinides.o: combineactinides.c header.h locations.h
	$(CC) $(CFLAGS) -c combineactinides.c

commandlineplotter.o: commandlineplotter.c header.h locations.h
	$(CC) $(CFLAGS) -c commandlineplotter.c

combinefissionyields.o: combinefissionyields.c header.h locations.h
	$(CC) $(CFLAGS) -c combinefissionyields.c

comparestr.o: comparestr.c header.h locations.h
	$(CC) $(CFLAGS) -c comparestr.c

complex.o: complex.c header.h locations.h
	$(CC) $(CFLAGS) -c complex.c

complexrea.o: complexrea.c header.h locations.h
	$(CC) $(CFLAGS) -c complexrea.c

contribsplit.o: contribsplit.c header.h locations.h
	$(CC) $(CFLAGS) -c contribsplit.c

coordexpans.o: coordexpans.c header.h locations.h
	$(CC) $(CFLAGS) -c coordexpans.c

coordtrans.o: coordtrans.c header.h locations.h
	$(CC) $(CFLAGS) -c coordtrans.c

countdynsrc.o: countdynsrc.c header.h locations.h
	$(CC) $(CFLAGS) -c countdynsrc.c

countsensparams.o: countsensparams.c header.h locations.h
	$(CC) $(CFLAGS) -c countsensparams.c

covmatrixfromblock.o: covmatrixfromblock.c header.h locations.h
	$(CC) $(CFLAGS) -c covmatrixfromblock.c

covmatrixfromsingle.o: covmatrixfromsingle.c header.h locations.h
	$(CC) $(CFLAGS) -c covmatrixfromsingle.c

createfinixifc.o: createfinixifc.c header.h locations.h
	$(CC) $(CFLAGS) -c createfinixifc.c

creategeometry.o: creategeometry.c header.h locations.h
	$(CC) $(CFLAGS) -c creategeometry.c

createmesh.o: createmesh.c header.h locations.h
	$(CC) $(CFLAGS) -c createmesh.c

createumshcells.o: createumshcells.c header.h locations.h
	$(CC) $(CFLAGS) -c createumshcells.c

createuniverse.o: createuniverse.c header.h locations.h
	$(CC) $(CFLAGS) -c createuniverse.c

cspline.o: cspline.c header.h
	$(CC) $(CFLAGS) -c cspline.c

cyldis.o: cyldis.c header.h
	$(CC) $(CFLAGS) -c cyldis.c

dataifcadens.o: dataifcadens.c header.h locations.h
	$(CC) $(CFLAGS) -c dataifcadens.c

dataifcxs.o: dataifcxs.c header.h locations.h
	$(CC) $(CFLAGS) -c dataifcxs.c

decaymeshprecdet.o: decaymeshprecdet.c header.h locations.h
	$(CC) $(CFLAGS) -c decaymeshprecdet.c

decaypointprecdet.o: decaypointprecdet.c header.h locations.h
	$(CC) $(CFLAGS) -c decaypointprecdet.c

decomposeelements.o: decomposeelements.c header.h locations.h element_data.h
	$(CC) $(CFLAGS) -c decomposeelements.c

defaultbradata.o: defaultbradata.c header.h locations.h
	$(CC) $(CFLAGS) -c defaultbradata.c

deinitsocket.o: deinitsocket.c header.h locations.h
	$(CC) $(CFLAGS) -c deinitsocket.c

densityfactor.o: densityfactor.c header.h locations.h
	$(CC) $(CFLAGS) -c densityfactor.c

depletionpolyfit.o: depletionpolyfit.c header.h locations.h
	$(CC) $(CFLAGS) -c depletionpolyfit.c

detbin.o: detbin.c header.h locations.h
	$(CC) $(CFLAGS) -c detbin.c

detectoroutput.o: detectoroutput.c header.h locations.h
	$(CC) $(CFLAGS) -c detectoroutput.c

deterministicleakage.o: deterministicleakage.c header.h locations.h
	$(CC) $(CFLAGS) -c deterministicleakage.c

detidx.o: detidx.c header.h locations.h
	$(CC) $(CFLAGS) -c detidx.c

detresponse.o: detresponse.c header.h locations.h
	$(CC) $(CFLAGS) -c detresponse.c

dfpos.o: dfpos.c header.h locations.h
	$(CC) $(CFLAGS) -c dfpos.c

dfsol.o: dfsol.c header.h locations.h
	$(CC) $(CFLAGS) -c dfsol.c

dfsolver.o: dfsolver.c header.h locations.h
	$(CC) $(CFLAGS) -c dfsolver.c

die.o: die.c header.h locations.h
	$(CC) $(CFLAGS) -c die.c

diffcoefed.o: diffcoefed.c header.h locations.h
	$(CC) $(CFLAGS) -c diffcoefed.c

disperse.o: disperse.c header.h locations.h
	$(CC) $(CFLAGS) -c disperse.c

disperse2.o: disperse2.c header.h locations.h
	$(CC) $(CFLAGS) -c disperse2.c

distributeddmajorants.o: distributeddmajorants.c header.h locations.h
	$(CC) $(CFLAGS) -c distributeddmajorants.c

distributeddmaxadens.o: distributeddmaxadens.c header.h locations.h
	$(CC) $(CFLAGS) -c distributeddmaxadens.c

distributeddlimbo.o: distributeddlimbo.c header.h locations.h
	$(CC) $(CFLAGS) -c distributeddlimbo.c

distributematerialdata.o: distributematerialdata.c header.h locations.h
	$(CC) $(CFLAGS) -c distributematerialdata.c

divideburnmat.o: divideburnmat.c header.h locations.h
	$(CC) $(CFLAGS) -c divideburnmat.c

dividemeshcell.o: dividemeshcell.c header.h locations.h
	$(CC) $(CFLAGS) -c dividemeshcell.c

dividepolyhedcell.o: dividepolyhedcell.c header.h locations.h
	$(CC) $(CFLAGS) -c dividepolyhedcell.c

dividepolyhedface.o: dividepolyhedface.c header.h locations.h
	$(CC) $(CFLAGS) -c dividepolyhedface.c

dividezone.o: dividezone.c header.h locations.h
	$(CC) $(CFLAGS) -c dividezone.c

divsubdistance.o: divsubdistance.c header.h locations.h
	$(CC) $(CFLAGS) -c divsubdistance.c

dopmicroxs.o: dopmicroxs.c header.h locations.h
	$(CC) $(CFLAGS) -c dopmicroxs.c

dopplerbroad.o: dopplerbroad.c header.h locations.h
	$(CC) $(CFLAGS) -c dopplerbroad.c

dtmajorant.o: dtmajorant.c header.h locations.h
	$(CC) $(CFLAGS) -c dtmajorant.c

duplicateitem.o: duplicateitem.c header.h locations.h
	$(CC) $(CFLAGS) -c duplicateitem.c

duplicateparticle.o: duplicateparticle.c header.h locations.h
	$(CC) $(CFLAGS) -c duplicateparticle.c

eblockfrombank.o: eblockfrombank.c header.h locations.h
	$(CC) $(CFLAGS) -c eblockfrombank.c

eblocktobank.o: eblocktobank.c header.h locations.h
	$(CC) $(CFLAGS) -c eblocktobank.c

elasticscattering.o: elasticscattering.c header.h locations.h
	$(CC) $(CFLAGS) -c elasticscattering.c

elbrsxsdata.o: elbrsxsdata.c header.h locations.h
	$(CC) $(CFLAGS) -c elbrsxsdata.c

element.o: element.c header.h natural_elements.h element_data.h locations.h
	$(CC) $(CFLAGS) -c element.c

elspdata.o: elspdata.c header.h locations.h
	$(CC) $(CFLAGS) -c elspdata.c

endfcolf.o: endfcolf.c header.h
	$(CC) $(CFLAGS) -c endfcolf.c

endfcoli.o: endfcoli.c header.h
	$(CC) $(CFLAGS) -c endfcoli.c

endfinterp.o: endfinterp.c header.h
	$(CC) $(CFLAGS) -c endfinterp.c

endfmat.o: endfmat.c header.h
	$(CC) $(CFLAGS) -c endfmat.c

endfmf.o: endfmf.c header.h
	$(CC) $(CFLAGS) -c endfmf.c

endfmt.o: endfmt.c header.h
	$(CC) $(CFLAGS) -c endfmt.c

endfnewline.o: endfnewline.c header.h
	$(CC) $(CFLAGS) -c endfnewline.c

eneimportance.o: eneimportance.c header.h locations.h input_params.h
	$(CC) $(CFLAGS) -c eneimportance.c

error.o: error.c header.h locations.h input_params.h
	$(CC) $(CFLAGS) -c error.c

eventfrombank.o: eventfrombank.c header.h locations.h
	$(CC) $(CFLAGS) -c eventfrombank.c

eventstosensitivity.o: eventstosensitivity.c header.h locations.h
	$(CC) $(CFLAGS) -c eventstosensitivity.c

eventtobank.o: eventtobank.c header.h locations.h
	$(CC) $(CFLAGS) -c eventtobank.c

estimateruntime.o: estimateruntime.c header.h locations.h
	$(CC) $(CFLAGS) -c estimateruntime.c

expandfe.o: expandfe.c header.h locations.h
	$(CC) $(CFLAGS) -c expandfe.c

expandprivatearrays.o: expandprivatearrays.c header.h locations.h
	$(CC) $(CFLAGS) -c expandprivatearrays.c

expodecfit.o: expodecfit.c header.h locations.h
	$(CC) $(CFLAGS) -c expodecfit.c

extractsenslabel.o: extractsenslabel.c header.h locations.h
	$(CC) $(CFLAGS) -c extractsenslabel.c

feifc.o: feifc.c header.h locations.h
	$(CC) $(CFLAGS) -c feifc.c

fetidx.o: fetidx.c header.h locations.h
	$(CC) $(CFLAGS) -c fetidx.c

fetfinalize.o: fetfinalize.c header.h locations.h
	$(CC) $(CFLAGS) -c fetfinalize.c

fillstlmesh.o: fillstlmesh.c header.h locations.h
	$(CC) $(CFLAGS) -c fillstlmesh.c

finalizeccprecdet.o: finalizeccprecdet.c header.h locations.h
	$(CC) $(CFLAGS) -c finalizeccprecdet.c

finalizempi.o: finalizempi.c header.h
	$(CC) $(CFLAGS) -c finalizempi.c

finalizeumsh.o: finalizeumsh.c header.h locations.h
	$(CC) $(CFLAGS) -c finalizeumsh.c

findregmeshifcindex.o: findregmeshifcindex.c header.h locations.h
	$(CC) $(CFLAGS) -c findregmeshifcindex.c

findsenseindex.o: findsenseindex.c header.h locations.h
	$(CC) $(CFLAGS) -c findsenseindex.c

findsensmuindex.o: findsensmuindex.c header.h locations.h
	$(CC) $(CFLAGS) -c findsensmuindex.c

findsenszaiindex.o: findsenszaiindex.c header.h locations.h
	$(CC) $(CFLAGS) -c findsenszaiindex.c

findsensindices.o: findsensindices.c header.h locations.h
	$(CC) $(CFLAGS) -c findsensindices.c

findtetcell.o: findtetcell.c header.h locations.h
	$(CC) $(CFLAGS) -c findtetcell.c

findinterfaceregions.o: findinterfaceregions.c header.h locations.h
	$(CC) $(CFLAGS) -c findinterfaceregions.c

findlatticeregion.o: findlatticeregion.c header.h locations.h
	$(CC) $(CFLAGS) -c findlatticeregion.c

findmaterialpointers.o: findmaterialpointers.c header.h locations.h
	$(CC) $(CFLAGS) -c findmaterialpointers.c

findnestregion.o: findnestregion.c header.h locations.h
	$(CC) $(CFLAGS) -c findnestregion.c

findnuclidedata.o: findnuclidedata.c header.h locations.h
	$(CC) $(CFLAGS) -c findnuclidedata.c

findrowindexes.o: findrowindexes.c header.h
	$(CC) $(CFLAGS) -c findrowindexes.c

findpbregion.o: findpbregion.c header.h locations.h
	$(CC) $(CFLAGS) -c findpbregion.c

findstlsolid.o: findstlsolid.c header.h locations.h
	$(CC) $(CFLAGS) -c findstlsolid.c

finduniversecell.o: finduniversecell.c header.h locations.h
	$(CC) $(CFLAGS) -c finduniversecell.c

findxslimits.o: findxslimits.c header.h
	$(CC) $(CFLAGS) -c findxslimits.c

finixptrfromfpe.o: finixptrfromfpe.c header.h locations.h
	$(CC) $(CFLAGS) -c finixptrfromfpe.c

firstitem.o: firstitem.c header.h locations.h
	$(CC) $(CFLAGS) -c firstitem.c

fisse.o: fisse.c header.h locations.h
	$(CC) $(CFLAGS) -c fisse.c

fissecomp.o: fissecomp.c header.h locations.h
	$(CC) $(CFLAGS) -c fissecomp.c

fissmtxindex.o: fissmtxindex.c header.h locations.h
	$(CC) $(CFLAGS) -c fissmtxindex.c

fissmtxoutput.o: fissmtxoutput.c header.h locations.h
	$(CC) $(CFLAGS) -c fissmtxoutput.c

fission.o: fission.c header.h locations.h
	$(CC) $(CFLAGS) -c fission.c

fixhexmesh.o: fixhexmesh.c header.h locations.h
	$(CC) $(CFLAGS) -c fixhexmesh.c

fixpolyhedmesh.o: fixpolyhedmesh.c header.h locations.h
	$(CC) $(CFLAGS) -c fixpolyhedmesh.c

flushbank.o: flushbank.c header.h locations.h
	$(CC) $(CFLAGS) -c flushbank.c

flushddparticles.o: flushddparticles.c header.h locations.h
	$(CC) $(CFLAGS) -c flushddparticles.c

flushprecsource.o: flushprecsource.c header.h locations.h
	$(CC) $(CFLAGS) -c flushprecsource.c

formtransmupaths.o: formtransmupaths.c header.h locations.h
	$(CC) $(CFLAGS) -c formtransmupaths.c

freeddcomms.o: freeddcomms.c header.h locations.h
	$(CC) $(CFLAGS) -c freeddcomms.c

freefinix.o: freefinix.c header.h locations.h
	$(CC) $(CFLAGS) -c freefinix.c

freemem.o: freemem.c header.h locations.h
	$(CC) $(CFLAGS) -c freemem.c

frombank.o: frombank.c header.h locations.h
	$(CC) $(CFLAGS) -c frombank.c

fromcommonque.o: fromcommonque.c header.h locations.h
	$(CC) $(CFLAGS) -c fromcommonque.c

fromsrc.o: fromsrc.c header.h locations.h
	$(CC) $(CFLAGS) -c fromsrc.c

fromstack.o: fromstack.c header.h locations.h
	$(CC) $(CFLAGS) -c fromstack.c

fromstore.o: fromstore.c header.h locations.h
	$(CC) $(CFLAGS) -c fromstore.c

fromtrkbank.o: fromtrkbank.c header.h locations.h
	$(CC) $(CFLAGS) -c fromtrkbank.c

fromque.o: fromque.c header.h locations.h
	$(CC) $(CFLAGS) -c fromque.c

gaussiansubst.o: gaussiansubst.c header.h
	$(CC) $(CFLAGS) -c gaussiansubst.c

geoimportance.o: geoimportance.c header.h locations.h
	$(CC) $(CFLAGS) -c geoimportance.c

geometryplotter.o: geometryplotter.c header.h locations.h
	$(CC) $(CFLAGS) -c geometryplotter.c

getbankedprecursors.o: getbankedprecursors.c header.h locations.h
	$(CC) $(CFLAGS) -c getbankedprecursors.c

getburnids.o: getburnids.c header.h locations.h
	$(CC) $(CFLAGS) -c getburnids.c

getifpancestordata.o: getifpancestordata.c header.h locations.h
	$(CC) $(CFLAGS) -c getifpancestordata.c

getrequiredcovmatrices.o: getrequiredcovmatrices.c header.h locations.h
	$(CC) $(CFLAGS) -c getrequiredcovmatrices.c

getsenseventtype.o: getsenseventtype.c header.h locations.h
	$(CC) $(CFLAGS) -c getsenseventtype.c

getimportantpts.o: getimportantpts.c header.h locations.h
	$(CC) $(CFLAGS) -c getimportantpts.c

getlatticeindexes.o: getlatticeindexes.c header.h locations.h
	$(CC) $(CFLAGS) -c getlatticeindexes.c

getlinenumber.o: getlinenumber.c header.h locations.h
	$(CC) $(CFLAGS) -c getlinenumber.c

getparams.o: getparams.c header.h input_params.h
	$(CC) $(CFLAGS) -c getparams.c

getprivatedata.o: getprivatedata.c header.h locations.h
	$(CC) $(CFLAGS) -c getprivatedata.c

getprivateres.o: getprivateres.c header.h locations.h
	$(CC) $(CFLAGS) -c getprivateres.c

getsignalfromsupervisor.o: getsignalfromsupervisor.c header.h locations.h
	$(CC) $(CFLAGS) -c getsignalfromsupervisor.c

gettemp.o: gettemp.c header.h locations.h
	$(CC) $(CFLAGS) -c gettemp.c

gettext.o: gettext.c header.h locations.h
	$(CC) $(CFLAGS) -c gettext.c

gridfactor.o: gridfactor.c header.h locations.h
	$(CC) $(CFLAGS) -c gridfactor.c

gridsearch.o: gridsearch.c header.h locations.h
	$(CC) $(CFLAGS) -c gridsearch.c

gsconfigdata.o: gsconfigdata.c  header.h locations.h
	$(CC) $(CFLAGS) -c gsconfigdata.c

hessfactorization.o: hessfactorization.c header.h locations.h
	$(CC) $(CFLAGS) -c hessfactorization.c

hexrotatecell.o: hexrotatecell.c header.h locations.h
	$(CC) $(CFLAGS) -c hexrotatecell.c

hexrotateface.o: hexrotateface.c header.h locations.h
	$(CC) $(CFLAGS) -c hexrotateface.c

homoflux.o: homoflux.c header.h locations.h
	$(CC) $(CFLAGS) -c homoflux.c

hismean.o: hismean.c header.h locations.h
	$(CC) $(CFLAGS) -c hismean.c

hisrelerr.o: hisrelerr.c header.h locations.h
	$(CC) $(CFLAGS) -c hisrelerr.c

hisval.o: hisval.c header.h locations.h
	$(CC) $(CFLAGS) -c hisval.c

icmidx.o: icmidx.c header.h locations.h
	$(CC) $(CFLAGS) -c icmidx.c

idxstr.o: idxstr.c header.h
	$(CC) $(CFLAGS) -c idxstr.c

ifcpoint.o: ifcpoint.c header.h locations.h
	$(CC) $(CFLAGS) -c ifcpoint.c

importancesolver.o: importancesolver.c header.h locations.h
	$(CC) $(CFLAGS) -c importancesolver.c

incell.o: incell.c header.h locations.h
	$(CC) $(CFLAGS) -c incell.c

inelasticscattering.o: inelasticscattering.c header.h locations.h
	$(CC) $(CFLAGS) -c inelasticscattering.c

initdata.o: initdata.c header.h locations.h
	$(CC) $(CFLAGS) -c initdata.c

initddcomms.o: initddcomms.c header.h locations.h
	$(CC) $(CFLAGS) -c initddcomms.c

initddrecv.o: initddrecv.c header.h locations.h
	$(CC) $(CFLAGS) -c initddrecv.c

inithistories.o: inithistories.c header.h locations.h
	$(CC) $(CFLAGS) -c inithistories.c

initialcritsrc.o: initialcritsrc.c header.h locations.h
	$(CC) $(CFLAGS) -c initialcritsrc.c

initmpi.o: initmpi.c header.h
	$(CC) $(CFLAGS) -c initmpi.c

initomp.o: initomp.c header.h locations.h
	$(CC) $(CFLAGS) -c initomp.c

initprecdet.o: initprecdet.c header.h locations.h
	$(CC) $(CFLAGS) -c initprecdet.c

initprecdetsource.o: initprecdetsource.c header.h locations.h
	$(CC) $(CFLAGS) -c initprecdetsource.c

initsignal.o: initsignal.c header.h locations.h
	$(CC) $(CFLAGS) -c initsignal.c

initsocket.o: initsocket.c header.h locations.h
	$(CC) $(CFLAGS) -c initsocket.c

inumshcell.o: inumshcell.c header.h locations.h
	$(CC) $(CFLAGS) -c inumshcell.c

insupercell.o: insupercell.c header.h locations.h
	$(CC) $(CFLAGS) -c insupercell.c

interpolatedata.o: interpolatedata.c header.h
	$(CC) $(CFLAGS) -c interpolatedata.c

interpolatefisse.o: interpolatefisse.c header.h
	$(CC) $(CFLAGS) -c interpolatefisse.c

interpolatemupdf.o: interpolatemupdf.c header.h
	$(CC) $(CFLAGS) -c interpolatemupdf.c

interpolatenubar.o: interpolatenubar.c header.h
	$(CC) $(CFLAGS) -c interpolatenubar.c

interpolatesab.o: interpolatesab.c header.h locations.h
	$(CC) $(CFLAGS) -c interpolatesab.c

intersectionlist.o: intersectionlist.c header.h locations.h
	$(CC) $(CFLAGS) -c intersectionlist.c

intetcell.o: intetcell.c header.h locations.h
	$(CC) $(CFLAGS) -c intetcell.c

invokebranch.o: invokebranch.c header.h locations.h
	$(CC) $(CFLAGS) -c invokebranch.c

involute.o: involute.c header.h
	$(CC) $(CFLAGS) -c involute.c

isotopefractions.o: isotopefractions.c header.h locations.h
	$(CC) $(CFLAGS) -c isotopefractions.c

isotozai.o: isotozai.c header.h element_data.h
	$(CC) $(CFLAGS) -c isotozai.c

isotropicdirection.o: isotropicdirection.c header.h
	$(CC) $(CFLAGS) -c isotropicdirection.c

iteratecc.o: iteratecc.c header.h locations.h
	$(CC) $(CFLAGS) -c iteratecc.c

iterateexternal.o: iterateexternal.c header.h locations.h
	$(CC) $(CFLAGS) -c iterateexternal.c

iteratefinix.o: iteratefinix.c header.h locations.h
	$(CC) $(CFLAGS) -c iteratefinix.c

iteratekeff.o: iteratekeff.c header.h locations.h
	$(CC) $(CFLAGS) -c iteratekeff.c

iternucxs.o: iternucxs.c header.h locations.h
	$(CC) $(CFLAGS) -c iternucxs.c

kleinnishina.o: kleinnishina.c header.h
	$(CC) $(CFLAGS) -c kleinnishina.c

lagrangeinterp.o: lagrangeinterp.c header.h
	$(CC) $(CFLAGS) -c lagrangeinterp.c

lastitem.o: lastitem.c header.h locations.h
	$(CC) $(CFLAGS) -c lastitem.c

latcellorig.o: latcellorig.c header.h locations.h
	$(CC) $(CFLAGS) -c latcellorig.c

leak.o: leak.c header.h locations.h
	$(CC) $(CFLAGS) -c leak.c

leakdet.o: leakdet.c header.h locations.h
	$(CC) $(CFLAGS) -c leakdet.c

levelscattering.o: levelscattering.c header.h locations.h
	$(CC) $(CFLAGS) -c levelscattering.c

lifolistsize.o: lifolistsize.c header.h locations.h
	$(CC) $(CFLAGS) -c lifolistsize.c

linkinterfacematerials.o: linkinterfacematerials.c header.h locations.h
	$(CC) $(CFLAGS) -c linkinterfacematerials.c

linkgcumaterials.o: linkgcumaterials.c header.h locations.h
	$(CC) $(CFLAGS) -c linkgcumaterials.c

linkreactions.o: linkreactions.c header.h locations.h
	$(CC) $(CFLAGS) -c linkreactions.c

linksabdata.o: linksabdata.c header.h locations.h
	$(CC) $(CFLAGS) -c linksabdata.c

listptr.o: listptr.c header.h locations.h
	$(CC) $(CFLAGS) -c listptr.c

ludecomposition.o: ludecomposition.c header.h locations.h
	$(CC) $(CFLAGS) -c ludecomposition.c

macrourescorr.o: macrourescorr.c header.h locations.h
	$(CC) $(CFLAGS) -c macrourescorr.c

macroxs.o: macroxs.c header.h locations.h
	$(CC) $(CFLAGS) -c macroxs.c

main.o: main.c header.h locations.h
	$(CC) $(CFLAGS) -c main.c

majorantxs.o: majorantxs.c header.h locations.h
	$(CC) $(CFLAGS) -c majorantxs.c

makearray.o: makearray.c header.h
	$(CC) $(CFLAGS) -c makearray.c

makeburnmatrix.o: makeburnmatrix.c header.h locations.h
	$(CC) $(CFLAGS) -c makeburnmatrix.c

makeburnmatrixmsr.o: makeburnmatrixmsr.c header.h locations.h
	$(CC) $(CFLAGS) -c makeburnmatrixmsr.c

makedepletionzones.o: makedepletionzones.c header.h locations.h
	$(CC) $(CFLAGS) -c makedepletionzones.c

makeenergygrid.o: makeenergygrid.c header.h locations.h
	$(CC) $(CFLAGS) -c makeenergygrid.c

makepalette.o: makepalette.c header.h locations.h
	$(CC) $(CFLAGS) -c makepalette.c

makering.o: makering.c header.h locations.h
	$(CC) $(CFLAGS) -c makering.c

materialburnup.o: materialburnup.c header.h locations.h
	$(CC) $(CFLAGS) -c materialburnup.c

materialtotals.o: materialtotals.c header.h locations.h
	$(CC) $(CFLAGS) -c materialtotals.c

materialvolumes.o: materialvolumes.c header.h locations.h
	$(CC) $(CFLAGS) -c materialvolumes.c

matproduct.o: matproduct.c header.h locations.h
	$(CC) $(CFLAGS) -c matproduct.c

matptr.o: matptr.c header.h locations.h
	$(CC) $(CFLAGS) -c matptr.c

matpos.o: matpos.c header.h locations.h
	$(CC) $(CFLAGS) -c matpos.c

matrixexponential.o: matrixexponential.c header.h locations.h
	$(CC) $(CFLAGS) -c matrixexponential.c

matlaboutput.o: matlaboutput.c header.h locations.h
	$(CC) $(CFLAGS) -c matlaboutput.c

maxsrcimp.o: maxsrcimp.c header.h locations.h
	$(CC) $(CFLAGS) -c maxsrcimp.c

maxsurfdimensions.o: maxsurfdimensions.c header.h locations.h
	$(CC) $(CFLAGS) -c maxsurfdimensions.c

maxwellenergy.o: maxwellenergy.c header.h
	$(CC) $(CFLAGS) -c maxwellenergy.c

mean.o: mean.c header.h locations.h
	$(CC) $(CFLAGS) -c mean.c

mem.o: mem.c header.h locations.h
	$(CC) $(CFLAGS) -c mem.c

memcount.o: memcount.c header.h locations.h
	$(CC) $(CFLAGS) -c memcount.c

meshcellbounds.o: meshcellbounds.c header.h locations.h
	$(CC) $(CFLAGS) -c meshcellbounds.c

meshcellconnectdiags.o: meshcellconnectdiags.c header.h locations.h
	$(CC) $(CFLAGS) -c meshcellconnectdiags.c

meshcellfromcgns.o: meshcellfromcgns.c header.h locations.h
	$(CC) $(CFLAGS) -c meshcellfromcgns.c

meshcellgetface.o: meshcellgetface.c header.h locations.h
	$(CC) $(CFLAGS) -c meshcellgetface.c

meshcellgetfacepos.o: meshcellgetfacepos.c header.h locations.h
	$(CC) $(CFLAGS) -c meshcellgetfacepos.c

meshcellindirection.o: meshcellindirection.c header.h locations.h
	$(CC) $(CFLAGS) -c meshcellindirection.c

meshcellrotatelists.o: meshcellrotatelists.c header.h locations.h
	$(CC) $(CFLAGS) -c meshcellrotatelists.c

meshcelltype.o: meshcelltype.c header.h locations.h
	$(CC) $(CFLAGS) -c meshcelltype.c

meshcellvol.o: meshcellvol.c header.h locations.h
	$(CC) $(CFLAGS) -c meshcellvol.c

meshindex.o: meshindex.c header.h locations.h
	$(CC) $(CFLAGS) -c meshindex.c

meshplotter.o: meshplotter.c header.h locations.h
	$(CC) $(CFLAGS) -c meshplotter.c

meshptr.o: meshptr.c header.h locations.h
	$(CC) $(CFLAGS) -c meshptr.c

meshtot.o: meshtot.c header.h locations.h
	$(CC) $(CFLAGS) -c meshtot.c

meshval.o: meshval.c header.h locations.h
	$(CC) $(CFLAGS) -c meshval.c

mgxs.o: mgxs.c header.h locations.h
	$(CC) $(CFLAGS) -c mgxs.c

microcalc.o: microcalc.c header.h locations.h
	$(CC) $(CFLAGS) -c microcalc.c

microdepoutput.o: microdepoutput.c header.h locations.h
	$(CC) $(CFLAGS) -c microdepoutput.c

micromajorantxs.o: micromajorantxs.c header.h locations.h
	$(CC) $(CFLAGS) -c micromajorantxs.c

microxs.o: microxs.c header.h locations.h
	$(CC) $(CFLAGS) -c microxs.c

minxs.o: minxs.c header.h locations.h
	$(CC) $(CFLAGS) -c minxs.c

moraoutput.o: moraoutput.c header.h locations.h
	$(CC) $(CFLAGS) -c moraoutput.c

movedt.o: movedt.c header.h locations.h
	$(CC) $(CFLAGS) -c movedt.c

moveifc.o: moveifc.c header.h locations.h
	$(CC) $(CFLAGS) -c moveifc.c

moveitemfirst.o: moveitemfirst.c header.h locations.h
	$(CC) $(CFLAGS) -c moveitemfirst.c

moveitemlast.o: moveitemlast.c header.h locations.h
	$(CC) $(CFLAGS) -c moveitemlast.c

moveitemright.o: moveitemright.c header.h locations.h
	$(CC) $(CFLAGS) -c moveitemright.c

movest.o: movest.c header.h locations.h
	$(CC) $(CFLAGS) -c movest.c

movestore.o: movestore.c header.h locations.h
	$(CC) $(CFLAGS) -c movestore.c

mpitransfer.o: mpitransfer.c header.h locations.h
	$(CC) $(CFLAGS) -c mpitransfer.c

msrrealist.o: msrrealist.c header.h locations.h
	$(CC) $(CFLAGS) -c msrrealist.c

myparallelmat.o: myparallelmat.c header.h locations.h
	$(CC) $(CFLAGS) -c myparallelmat.c

nearestboundary.o: nearestboundary.c header.h locations.h
	$(CC) $(CFLAGS) -c nearestboundary.c

nearestpbsurf.o: nearestpbsurf.c header.h locations.h
	$(CC) $(CFLAGS) -c nearestpbsurf.c

nearestmeshboundary.o: nearestmeshboundary.c header.h locations.h
	$(CC) $(CFLAGS) -c nearestmeshboundary.c

neareststlsurf.o: neareststlsurf.c header.h locations.h
	$(CC) $(CFLAGS) -c neareststlsurf.c

nearestumshsurf.o: nearestumshsurf.c header.h locations.h
	$(CC) $(CFLAGS) -c nearestumshsurf.c

nestvolumes.o: nestvolumes.c header.h locations.h
	$(CC) $(CFLAGS) -c nestvolumes.c

newitem.o: newitem.c header.h locations.h
	$(CC) $(CFLAGS) -c newitem.c

newlifoitem.o: newlifoitem.c header.h locations.h
	$(CC) $(CFLAGS) -c newlifoitem.c

newrealist.o: newrealist.c header.h locations.h
	$(CC) $(CFLAGS) -c newrealist.c

newstat.o: newstat.c header.h locations.h
	$(CC) $(CFLAGS) -c newstat.c

nextitem.o: nextitem.c header.h locations.h
	$(CC) $(CFLAGS) -c nextitem.c

nextreaction.o: nextreaction.c header.h locations.h
	$(CC) $(CFLAGS) -c nextreaction.c

nextword.o: nextword.c header.h
	$(CC) $(CFLAGS) -c nextword.c

nfkerma.o: nfkerma.c header.h locations.h
	$(CC) $(CFLAGS) -c nfkerma.c

nocorrector.o: nocorrector.c header.h locations.h
	$(CC) $(CFLAGS) -c nocorrector.c

normcoef.o: normcoef.c header.h locations.h
	$(CC) $(CFLAGS) -c normcoef.c

normalizecompositions.o: normalizecompositions.c header.h locations.h
	$(CC) $(CFLAGS) -c normalizecompositions.c

normalizecritsrc.o: normalizecritsrc.c header.h locations.h
	$(CC) $(CFLAGS) -c normalizecritsrc.c

normalizedynsrc.o: normalizedynsrc.c header.h locations.h
	$(CC) $(CFLAGS) -c normalizedynsrc.c

normalizeimp.o: normalizeimp.c header.h locations.h
	$(CC) $(CFLAGS) -c normalizeimp.c

normalizeprecdet.o: normalizeprecdet.c header.h locations.h
	$(CC) $(CFLAGS) -c normalizeprecdet.c

note.o: note.c header.h
	$(CC) $(CFLAGS) -c note.c

nubar.o: nubar.c header.h locations.h
	$(CC) $(CFLAGS) -c nubar.c

numericgauss.o: numericgauss.c header.h
	$(CC) $(CFLAGS) -c numericgauss.c

numericstr.o: numericstr.c header.h
	$(CC) $(CFLAGS) -c numericstr.c

nxn.o: nxn.c header.h locations.h
	$(CC) $(CFLAGS) -c nxn.c

ompresetcomp.o: ompresetcomp.c header.h locations.h
	$(CC) $(CFLAGS) -c ompresetcomp.c

ompsetcomp.o: ompsetcomp.c header.h locations.h
	$(CC) $(CFLAGS) -c ompsetcomp.c

omptestcomp.o: omptestcomp.c header.h locations.h
	$(CC) $(CFLAGS) -c omptestcomp.c

opendatafile.o: opendatafile.c header.h
	$(CC) $(CFLAGS) -c opendatafile.c

otfburntta.o: otfburntta.c header.h locations.h
	$(CC) $(CFLAGS) -c otfburntta.c

otfburnxs.o: otfburnxs.c header.h locations.h
	$(CC) $(CFLAGS) -c otfburnxs.c

otfsabxs.o: otfsabxs.c header.h locations.h
	$(CC) $(CFLAGS) -c otfsabxs.c

otfsabscattering.o: otfsabscattering.c header.h locations.h
	$(CC) $(CFLAGS) -c otfsabscattering.c

overrideids.o: overrideids.c header.h locations.h
	$(CC) $(CFLAGS) -c overrideids.c

pairproduction.o: pairproduction.c header.h locations.h
	$(CC) $(CFLAGS) -c pairproduction.c

parlett.o: parlett.c header.h locations.h
	$(CC) $(CFLAGS) -c parlett.c

parsecommandline.o: parsecommandline.c header.h locations.h
	$(CC) $(CFLAGS) -c parsecommandline.c

particlesfromstore.o: particlesfromstore.c header.h locations.h
	$(CC) $(CFLAGS) -c particlesfromstore.c

pause.o: pause.c header.h
	$(CC) $(CFLAGS) -c pause.c

phidis.o: phidis.c header.h
	$(CC) $(CFLAGS) -c phidis.c

photoelectric.o: photoelectric.c header.h locations.h
	$(CC) $(CFLAGS) -c photoelectric.c

photonmacroxs.o: photonmacroxs.c header.h locations.h
	$(CC) $(CFLAGS) -c photonmacroxs.c

photonmicroxs.o: photonmicroxs.c header.h locations.h
	$(CC) $(CFLAGS) -c photonmicroxs.c

photonprod.o: photonprod.c header.h locations.h
	$(CC) $(CFLAGS) -c photonprod.c

plotimage.o: plotimage.c header.h locations.h
	$(CC) $(CFLAGS) -c plotimage.c

plottracks.o: plottracks.c header.h locations.h
	$(CC) $(CFLAGS) -c plottracks.c

poisoneq.o: poisoneq.c header.h locations.h
	$(CC) $(CFLAGS) -c poisoneq.c

poisonxs.o: poisonxs.c header.h locations.h
	$(CC) $(CFLAGS) -c poisonxs.c

polarangle.o: polarangle.c header.h locations.h
	$(CC) $(CFLAGS) -c polarangle.c

polynomiallegendre.o: polynomiallegendre.c header.h
	$(CC) $(CFLAGS) -c polynomiallegendre.c

polynomialzernike.o: polynomialzernike.c header.h
	$(CC) $(CFLAGS) -c polynomialzernike.c

polypinf.o: polypinf.c header.h
	$(CC) $(CFLAGS) -c polypinf.c

polysameface.o: polysameface.c header.h
	$(CC) $(CFLAGS) -c polysameface.c

posannihilation.o: posannihilation.c header.h locations.h
	$(CC) $(CFLAGS) -c posannihilation.c

postddrecv.o: postddrecv.c header.h locations.h
	$(CC) $(CFLAGS) -c postddrecv.c

postddsend.o: postddsend.c header.h locations.h
	$(CC) $(CFLAGS) -c postddsend.c

potcorr.o: potcorr.c header.h locations.h
	$(CC) $(CFLAGS) -c potcorr.c

preallocmem.o: preallocmem.c header.h locations.h
	$(CC) $(CFLAGS) -c preallocmem.c

precdet.o: precdet.c header.h locations.h
	$(CC) $(CFLAGS) -c precdet.c

precursorpopcontrol.o: precursorpopcontrol.c header.h locations.h
	$(CC) $(CFLAGS) -c precursorpopcontrol.c

precursorstostore.o: precursorstostore.c header.h locations.h
	$(CC) $(CFLAGS) -c precursorstostore.c

preparecciter.o: preparecciter.c header.h locations.h
	$(CC) $(CFLAGS) -c preparecciter.c

preparecellsearchmesh.o: preparecellsearchmesh.c header.h locations.h
	$(CC) $(CFLAGS) -c preparecellsearchmesh.c

preparetransportcycle.o: preparetransportcycle.c header.h locations.h
	$(CC) $(CFLAGS) -c preparetransportcycle.c

presort.o: presort.c header.h locations.h
	$(CC) $(CFLAGS) -c presort.c

pretrans.o: pretrans.c header.h locations.h
	$(CC) $(CFLAGS) -c pretrans.c

previtem.o: previtem.c header.h locations.h
	$(CC) $(CFLAGS) -c previtem.c

printcellstats.o: printcellstats.c header.h locations.h
	$(CC) $(CFLAGS) -c printcellstats.c

printcoevals.o: printcoevals.c header.h locations.h
	$(CC) $(CFLAGS) -c printcoevals.c

printcompositions.o: printcompositions.c header.h locations.h
	$(CC) $(CFLAGS) -c printcompositions.c

printcoredistr.o: printcoredistr.c header.h locations.h
	$(CC) $(CFLAGS) -c printcoredistr.c

printcycleoutput.o: printcycleoutput.c header.h locations.h
	$(CC) $(CFLAGS) -c printcycleoutput.c

printdepmatrix.o: printdepmatrix.c header.h locations.h
	$(CC) $(CFLAGS) -c printdepmatrix.c

printdepoutput.o: printdepoutput.c header.h locations.h
	$(CC) $(CFLAGS) -c printdepoutput.c

printdepvals.o: printdepvals.c header.h locations.h
	$(CC) $(CFLAGS) -c printdepvals.c

printfinix.o: printfinix.c header.h locations.h
	$(CC) $(CFLAGS) -c printfinix.c

printgammaspectra.o: printgammaspectra.c header.h locations.h
	$(CC) $(CFLAGS) -c printgammaspectra.c

printgeometrydata.o: printgeometrydata.c header.h locations.h
	$(CC) $(CFLAGS) -c printgeometrydata.c

printhistoryoutput.o: printhistoryoutput.c header.h locations.h
	$(CC) $(CFLAGS) -c printhistoryoutput.c

printinterfaceoutput.o: printinterfaceoutput.c header.h locations.h
	$(CC) $(CFLAGS) -c printinterfaceoutput.c

printmaterialdata.o: printmaterialdata.c header.h locations.h
	$(CC) $(CFLAGS) -c printmaterialdata.c

printmeshcell.o: printmeshcell.c header.h locations.h
	$(CC) $(CFLAGS) -c printmeshcell.c

printmixtures.o: printmixtures.c header.h locations.h
	$(CC) $(CFLAGS) -c printmixtures.c

printmvar.o: printmvar.c header.h locations.h
	$(CC) $(CFLAGS) -c printmvar.c

printnuclidedata.o: printnuclidedata.c header.h locations.h
	$(CC) $(CFLAGS) -c printnuclidedata.c

printpbdata.o: printpbdata.c header.h locations.h
	$(CC) $(CFLAGS) -c printpbdata.c

printprecdet.o: printprecdet.c header.h locations.h
	$(CC) $(CFLAGS) -c printprecdet.c

printprogress.o: printprogress.c header.h locations.h
	$(CC) $(CFLAGS) -c printprogress.c

printsrcimp.o: printsrcimp.c header.h locations.h
	$(CC) $(CFLAGS) -c printsrcimp.c

printreactionlists.o: printreactionlists.c header.h locations.h
	$(CC) $(CFLAGS) -c printreactionlists.c

printtitle.o: printtitle.c header.h locations.h
	$(CC) $(CFLAGS) -c printtitle.c

printtmsdiagnostics.o: printtmsdiagnostics.c header.h locations.h
	$(CC) $(CFLAGS) -c printtmsdiagnostics.c

printvalues.o: printvalues.c header.h locations.h
	$(CC) $(CFLAGS) -c printvalues.c

probevrmesh.o: probevrmesh.c header.h locations.h
	$(CC) $(CFLAGS) -c probevrmesh.c

processbc.o: processbc.c header.h locations.h
	$(CC) $(CFLAGS) -c processbc.c

processbuinterface.o: processbuinterface.c header.h locations.h
	$(CC) $(CFLAGS) -c processbuinterface.c

processburnmat.o: processburnmat.c header.h locations.h
	$(CC) $(CFLAGS) -c processburnmat.c

processburnupegroups.o: processburnupegroups.c header.h locations.h
	$(CC) $(CFLAGS) -c processburnupegroups.c

processcells.o: processcells.c header.h locations.h
	$(CC) $(CFLAGS) -c processcells.c

processcovariancedata.o: processcovariancedata.c header.h locations.h
	$(CC) $(CFLAGS) -c processcovariancedata.c

processcellmesh.o: processcellmesh.c header.h locations.h
	$(CC) $(CFLAGS) -c processcellmesh.c

processcpd.o: processcpd.c header.h locations.h
	$(CC) $(CFLAGS) -c processcpd.c

processcomplementcells.o: processcomplementcells.c header.h locations.h
	$(CC) $(CFLAGS) -c processcomplementcells.c

processdatainterfaces.o: processdatainterfaces.c header.h locations.h
	$(CC) $(CFLAGS) -c processdatainterfaces.c

processcompton.o: processcompton.c header.h locations.h
	$(CC) $(CFLAGS) -c processcompton.c

processdecaydata.o: processdecaydata.c header.h locations.h radiotox.h
	$(CC) $(CFLAGS) -c processdecaydata.c

processdecaysrc.o: processdecaysrc.c header.h locations.h
	$(CC) $(CFLAGS) -c processdecaysrc.c

processdephis.o: processdephis.c header.h locations.h
	$(CC) $(CFLAGS) -c processdephis.c

processdetectors.o: processdetectors.c header.h locations.h photon_attenuation.h
	$(CC) $(CFLAGS) -c processdetectors.c

processdivisors.o: processdivisors.c header.h locations.h
	$(CC) $(CFLAGS) -c processdivisors.c

processdix.o: processdix.c header.h locations.h
	$(CC) $(CFLAGS) -c processdix.c

processentropy.o: processentropy.c header.h locations.h
	$(CC) $(CFLAGS) -c processentropy.c

processedistributions.o: processedistributions.c header.h locations.h
	$(CC) $(CFLAGS) -c processedistributions.c

processelectrons.o: processelectrons.c header.h locations.h
	$(CC) $(CFLAGS) -c processelectrons.c

processevents.o: processevents.c header.h locations.h
	$(CC) $(CFLAGS) -c processevents.c

processfinix.o: processfinix.c header.h locations.h
	$(CC) $(CFLAGS) -c processfinix.c

processfissmtx.o: processfissmtx.c header.h locations.h
	$(CC) $(CFLAGS) -c processfissmtx.c

processfissionyields.o: processfissionyields.c header.h locations.h
	$(CC) $(CFLAGS) -c processfissionyields.c

processgc.o: processgc.c header.h locations.h
	$(CC) $(CFLAGS) -c processgc.c

processifcfb.o: processifcfb.c header.h locations.h
	$(CC) $(CFLAGS) -c processifcfb.c

processifcfetmaterial.o: processifcfetmaterial.c header.h locations.h
	$(CC) $(CFLAGS) -c processifcfetmaterial.c

processifcfunc.o: processifcfunc.c header.h locations.h
	$(CC) $(CFLAGS) -c processifcfunc.c

processifcptavg.o: processifcptavg.c header.h locations.h
	$(CC) $(CFLAGS) -c processifcptavg.c

processifcregmesh.o: processifcregmesh.c header.h locations.h
	$(CC) $(CFLAGS) -c processifcregmesh.c

processifctetmesh.o: processifctetmesh.c header.h locations.h
	$(CC) $(CFLAGS) -c processifctetmesh.c

processicm.o: processicm.c header.h locations.h
	$(CC) $(CFLAGS) -c processicm.c

processinterface.o: processinterface.c header.h locations.h
	$(CC) $(CFLAGS) -c processinterface.c

processinelaxs.o: processinelaxs.c header.h locations.h
	$(CC) $(CFLAGS) -c processinelaxs.c

processinventory.o: processinventory.c header.h locations.h
	$(CC) $(CFLAGS) -c processinventory.c

processiternucs.o: processiternucs.c header.h locations.h
	$(CC) $(CFLAGS) -c processiternucs.c

processlattices.o: processlattices.c header.h locations.h
	$(CC) $(CFLAGS) -c processlattices.c

processmaterials.o: processmaterials.c header.h locations.h
	$(CC) $(CFLAGS) -c processmaterials.c

processmeshplots.o: processmeshplots.c header.h locations.h
	$(CC) $(CFLAGS) -c processmeshplots.c

processmixture.o: processmixture.c header.h locations.h
	$(CC) $(CFLAGS) -c processmixture.c

processmsr.o: processmsr.c header.h locations.h
	$(CC) $(CFLAGS) -c processmsr.c

processmudistributions.o: processmudistributions.c header.h locations.h
	$(CC) $(CFLAGS) -c processmudistributions.c

processnests.o: processnests.c header.h locations.h
	$(CC) $(CFLAGS) -c processnests.c

processnubardata.o: processnubardata.c header.h locations.h
	$(CC) $(CFLAGS) -c processnubardata.c

processnuclides.o: processnuclides.c header.h locations.h
	$(CC) $(CFLAGS) -c processnuclides.c

processotfburn.o: processotfburn.c header.h locations.h
	$(CC) $(CFLAGS) -c processotfburn.c

processpairproduction.o: processpairproduction.c header.h locations.h
	$(CC) $(CFLAGS) -c processpairproduction.c

processpbgeometry.o: processpbgeometry.c header.h locations.h
	$(CC) $(CFLAGS) -c processpbgeometry.c

processphotoelectric.o: processphotoelectric.c header.h locations.h
	$(CC) $(CFLAGS) -c processphotoelectric.c

processphotonatt.o: processphotonatt.c header.h locations.h photon_attenuation.h
	$(CC) $(CFLAGS) -c processphotonatt.c

processphotonecut.o: processphotonecut.c header.h locations.h
	$(CC) $(CFLAGS) -c processphotonecut.c

processphotonprod.o: processphotonprod.c header.h locations.h
	$(CC) $(CFLAGS) -c processphotonprod.c

processphotonrea.o: processphotonrea.c header.h locations.h
	$(CC) $(CFLAGS) -c processphotonrea.c

processpoisons.o: processpoisons.c header.h locations.h
	$(CC) $(CFLAGS) -c processpoisons.c

processprecdet.o: processprecdet.c header.h locations.h
	$(CC) $(CFLAGS) -c processprecdet.c

processrayleigh.o: processrayleigh.c header.h locations.h
	$(CC) $(CFLAGS) -c processrayleigh.c

processreactionlists.o: processreactionlists.c header.h locations.h
	$(CC) $(CFLAGS) -c processreactionlists.c

processrelaxation.o : processrelaxation.c header.h locations.h
	$(CC) $(CFLAGS) -c processrelaxation.c

processreprocessors.o: processreprocessors.c header.h locations.h
	$(CC) $(CFLAGS) -c processreprocessors.c

processrmx.o: processrmx.c header.h locations.h
	$(CC) $(CFLAGS) -c processrmx.c

processsenseblocks.o: processsenseblocks.c header.h locations.h
	$(CC) $(CFLAGS) -c processsenseblocks.c

processsensitivities.o: processsensitivities.c header.h locations.h
	$(CC) $(CFLAGS) -c processsensitivities.c

processsensstats.o: processsensstats.c header.h locations.h
	$(CC) $(CFLAGS) -c processsensstats.c

processsingleinterface.o: processsingleinterface.c header.h locations.h
	$(CC) $(CFLAGS) -c processsingleinterface.c

processsources.o: processsources.c header.h locations.h
	$(CC) $(CFLAGS) -c processsources.c

processstats.o: processstats.c header.h locations.h
	$(CC) $(CFLAGS) -c processstats.c

processstlgeometry.o: processstlgeometry.c header.h locations.h
	$(CC) $(CFLAGS) -c processstlgeometry.c

processsurfaces.o: processsurfaces.c header.h locations.h
	$(CC) $(CFLAGS) -c processsurfaces.c

processsymmetries.o: processsymmetries.c header.h locations.h
	$(CC) $(CFLAGS) -c processsymmetries.c

processtmpdata.o: processtmpdata.c header.h locations.h
	$(CC) $(CFLAGS) -c processtmpdata.c

processtimebins.o: processtimebins.c header.h locations.h
	$(CC) $(CFLAGS) -c processtimebins.c

processtransformations.o: processtransformations.c header.h locations.h
	$(CC) $(CFLAGS) -c processtransformations.c

processtranspcorr.o: processtranspcorr.c header.h locations.h
	$(CC) $(CFLAGS) -c processtranspcorr.c

processttb.o: processttb.c header.h locations.h
	$(CC) $(CFLAGS) -c processttb.c

processumshgeometry.o: processumshgeometry.c header.h locations.h
	$(CC) $(CFLAGS) -c processumshgeometry.c

processuresdata.o: processuresdata.c header.h locations.h
	$(CC) $(CFLAGS) -c processuresdata.c

processuseregrids.o: processuseregrids.c header.h locations.h
	$(CC) $(CFLAGS) -c processuseregrids.c

processvr.o: processvr.c header.h locations.h
	$(CC) $(CFLAGS) -c processvr.c

processxsdata.o: processxsdata.c header.h locations.h
	$(CC) $(CFLAGS) -c processxsdata.c

pulsedet.o: pulsedet.c header.h locations.h
	$(CC) $(CFLAGS) -c pulsedet.c

putprivatedata.o: putprivatedata.c header.h locations.h
	$(CC) $(CFLAGS) -c putprivatedata.c

putcompositions.o: putcompositions.c header.h locations.h
	$(CC) $(CFLAGS) -c putcompositions.c

putmeshidx.o: putmeshidx.c header.h locations.h
	$(CC) $(CFLAGS) -c putmeshidx.c

putotfburnconc.o: putotfburnconc.c header.h locations.h
	$(CC) $(CFLAGS) -c putotfburnconc.c

putplotcolors.o: putplotcolors.c header.h locations.h
	$(CC) $(CFLAGS) -c putplotcolors.c

putpoisonconc.o: putpoisonconc.c header.h locations.h
	$(CC) $(CFLAGS) -c putpoisonconc.c

puttext.o: puttext.c header.h locations.h
	$(CC) $(CFLAGS) -c puttext.c

pythonplotter.o: pythonplotter.c header.h locations.h
	$(CC) $(CFLAGS) -c pythonplotter.c

qrfactorization.o: qrfactorization.c header.h locations.h
	$(CC) $(CFLAGS) -c qrfactorization.c

radgammasrc.o: radgammasrc.c header.h locations.h
	$(CC) $(CFLAGS) -c radgammasrc.c

rand64.o: rand64.c header.h
	$(CC) $(CFLAGS) -c rand64.c

randf.o: randf.c header.h locations.h
	$(CC) $(CFLAGS) -c randf.c

randnorm.o: randnorm.c header.h 
	$(CC) $(CFLAGS) -c randnorm.c

rayleighscattering.o: rayleighscattering.c header.h locations.h
	$(CC) $(CFLAGS) -c rayleighscattering.c

reactioncount.o: reactioncount.c header.h locations.h
	$(CC) $(CFLAGS) -c reactioncount.c

reactioncutoff.o: reactioncutoff.c header.h locations.h
	$(CC) $(CFLAGS) -c reactioncutoff.c

reactionmt.o: reactionmt.c header.h
	$(CC) $(CFLAGS) -c reactionmt.c

reactiontargetzai.o: reactiontargetzai.c header.h locations.h
	$(CC) $(CFLAGS) -c reactiontargetzai.c

readacefile.o: readacefile.c header.h locations.h
	$(CC) $(CFLAGS) -c readacefile.c

readdecayfile.o: readdecayfile.c header.h locations.h
	$(CC) $(CFLAGS) -c readdecayfile.c

readcoverxfile.o: readcoverxfile.c header.h locations.h
	$(CC) $(CFLAGS) -c readcoverxfile.c

readdatainterfaces.o: readdatainterfaces.c header.h locations.h
	$(CC) $(CFLAGS) -c readdatainterfaces.c

readdirectoryfile.o: readdirectoryfile.c header.h locations.h
	$(CC) $(CFLAGS) -c readdirectoryfile.c

readfinixifc.o: readfinixifc.c header.h locations.h
	$(CC) $(CFLAGS) -c readfinixifc.c

readfissionyields.o: readfissionyields.c header.h locations.h
	$(CC) $(CFLAGS) -c readfissionyields.c

readbrafile.o: readbrafile.c header.h locations.h
	$(CC) $(CFLAGS) -c readbrafile.c

readifcbins.o: readifcbins.c header.h locations.h
	$(CC) $(CFLAGS) -c readifcbins.c

readifcfb.o: readifcfb.c header.h locations.h
	$(CC) $(CFLAGS) -c readifcfb.c

readifcfblims.o: readifcfblims.c header.h locations.h
	$(CC) $(CFLAGS) -c readifcfblims.c

readifcfe.o: readifcfe.c header.h locations.h
	$(CC) $(CFLAGS) -c readifcfe.c

readifcfunc.o: readifcfunc.c header.h locations.h
	$(CC) $(CFLAGS) -c readifcfunc.c

readifcofmesh.o: readifcofmesh.c header.h locations.h
	$(CC) $(CFLAGS) -c readifcofmesh.c

readifcptavg.o: readifcptavg.c header.h locations.h
	$(CC) $(CFLAGS) -c readifcptavg.c

readifcregmesh.o: readifcregmesh.c header.h locations.h
	$(CC) $(CFLAGS) -c readifcregmesh.c

readifctetmesh.o: readifctetmesh.c header.h locations.h
	$(CC) $(CFLAGS) -c readifctetmesh.c

readinfix.o: readinfix.c header.h locations.h
	$(CC) $(CFLAGS) -c readinfix.c

readinput.o: readinput.c header.h locations.h surface_types.h
	$(CC) $(CFLAGS) -c readinput.c

readinterface.o: readinterface.c header.h locations.h
	$(CC) $(CFLAGS) -c readinterface.c

readmesh.o: readmesh.c header.h locations.h
	$(CC) $(CFLAGS) -c readmesh.c

readmeshptr.o: readmeshptr.c header.h locations.h
	$(CC) $(CFLAGS) -c readmeshptr.c

readpendfdata.o: readpendfdata.c header.h locations.h
	$(CC) $(CFLAGS) -c readpendfdata.c

readofpatches.o: readofpatches.c header.h locations.h
	$(CC) $(CFLAGS) -c readofpatches.c

readofdata.o: readofdata.c header.h  locations.h
	$(CC) $(CFLAGS) -c readofdata.c

readofdensities.o: readofdensities.c header.h  locations.h
	$(CC) $(CFLAGS) -c readofdensities.c

readoffaces.o: readoffaces.c header.h  locations.h
	$(CC) $(CFLAGS) -c readoffaces.c

readofheader.o: readofheader.c header.h locations.h
	$(CC) $(CFLAGS) -c readofheader.c

readofmapping.o: readofmapping.c header.h locations.h
	$(CC) $(CFLAGS) -c readofmapping.c

readofmaterials.o: readofmaterials.c header.h locations.h
	$(CC) $(CFLAGS) -c readofmaterials.c

readofmesh.o: readofmesh.c header.h locations.h
	$(CC) $(CFLAGS) -c readofmesh.c

readofownnbr.o: readofownnbr.c header.h locations.h
	$(CC) $(CFLAGS) -c readofownnbr.c

readofpoints.o: readofpoints.c header.h locations.h
	$(CC) $(CFLAGS) -c readofpoints.c

readoftemperatures.o: readoftemperatures.c header.h locations.h
	$(CC) $(CFLAGS) -c readoftemperatures.c

readpbgeometry.o: readpbgeometry.c header.h locations.h
	$(CC) $(CFLAGS) -c readpbgeometry.c

readplasmasrc.o: readplasmasrc.c header.h locations.h
	$(CC) $(CFLAGS) -c readplasmasrc.c

readphotondata.o: readphotondata.c header.h locations.h
	$(CC) $(CFLAGS) -c readphotondata.c

readrestartfile.o: readrestartfile.c header.h locations.h
	$(CC) $(CFLAGS) -c readrestartfile.c

readsourcefile.o: readsourcefile.c header.h locations.h
	$(CC) $(CFLAGS) -c readsourcefile.c

readstlgeometry.o: readstlgeometry.c header.h locations.h
	$(CC) $(CFLAGS) -c readstlgeometry.c

readstlmesh.o: readstlmesh.c header.h locations.h
	$(CC) $(CFLAGS) -c readstlmesh.c

readtextfile.o: readtextfile.c header.h locations.h 
	$(CC) $(CFLAGS) -c readtextfile.c

readumshgeometry.o: readumshgeometry.c header.h locations.h
	$(CC) $(CFLAGS) -c readumshgeometry.c

readwwdmesh.o: readwwdmesh.c header.h locations.h
	$(CC) $(CFLAGS) -c readwwdmesh.c

reallocddsend.o: reallocddsend.c header.h locations.h
	$(CC) $(CFLAGS) -c reallocddsend.c

reallocmem.o: reallocmem.c header.h locations.h
	$(CC) $(CFLAGS) -c reallocmem.c

reamulti.o: reamulti.c header.h locations.h
	$(CC) $(CFLAGS) -c reamulti.c

receiveifcinputdata.o: receiveifcinputdata.c header.h locations.h
	$(CC) $(CFLAGS) -c receiveifcinputdata.c

receiveddparticlecounters.o: receiveddparticlecounters.c header.h locations.h
	$(CC) $(CFLAGS) -c receiveddparticlecounters.c

receiveddparticles.o: receiveddparticles.c header.h locations.h
	$(CC) $(CFLAGS) -c receiveddparticles.c

recoildet.o: recoildet.c header.h locations.h
	$(CC) $(CFLAGS) -c recoildet.c

recyclermxdata.o: recyclermxdata.c header.h locations.h
	$(CC) $(CFLAGS) -c recyclermxdata.c

redistributeques.o: redistributeques.c header.h locations.h
	$(CC) $(CFLAGS) -c redistributeques.c

redistributestacks.o: redistributestacks.c header.h locations.h
	$(CC) $(CFLAGS) -c redistributestacks.c

reducebuffer.o: reducebuffer.c header.h locations.h
	$(CC) $(CFLAGS) -c reducebuffer.c

reduceprivateres.o: reduceprivateres.c header.h locations.h
	$(CC) $(CFLAGS) -c reduceprivateres.c

refreshinventory.o: refreshinventory.c header.h locations.h
	$(CC) $(CFLAGS) -c refreshinventory.c

reinitrng.o: reinitrng.c header.h locations.h
	$(CC) $(CFLAGS) -c reinitrng.c

relaxccresults.o: relaxccresults.c header.h locations.h
	$(CC) $(CFLAGS) -c relaxccresults.c

relaxinterfacepower.o: relaxinterfacepower.c header.h locations.h
	$(CC) $(CFLAGS) -c relaxinterfacepower.c

relaxtransmuxs.o: relaxtransmuxs.c header.h locations.h
	$(CC) $(CFLAGS) -c relaxtransmuxs.c

relerr.o: relerr.c header.h locations.h
	$(CC) $(CFLAGS) -c relerr.c

removeflaggeditems.o: removeflaggeditems.c header.h locations.h
	$(CC) $(CFLAGS) -c removeflaggeditems.c

removeitem.o: removeitem.c header.h locations.h
	$(CC) $(CFLAGS) -c removeitem.c

removestorefiles.o: removestorefiles.c header.h locations.h
	$(CC) $(CFLAGS) -c removestorefiles.c

removevoidcells.o: removevoidcells.c header.h locations.h
	$(CC) $(CFLAGS) -c removevoidcells.c

reopenlist.o: reopenlist.c header.h locations.h
	$(CC) $(CFLAGS) -c reopenlist.c

replaceitem.o: replaceitem.c header.h locations.h
	$(CC) $(CFLAGS) -c replaceitem.c

replacephotondata.o: replacephotondata.c header.h locations.h
	$(CC) $(CFLAGS) -c replacephotondata.c

reprocess.o: reprocess.c header.h locations.h
	$(CC) $(CFLAGS) -c reprocess.c

resetddcomms.o: resetddcomms.c header.h locations.h
	$(CC) $(CFLAGS) -c resetddcomms.c

resetddsend.o: resetddsend.c header.h locations.h
	$(CC) $(CFLAGS) -c resetddsend.c

resetoption.o: resetoption.c header.h locations.h
	$(CC) $(CFLAGS) -c resetoption.c

resetotfburn.o: resetotfburn.c header.h locations.h
	$(CC) $(CFLAGS) -c resetotfburn.c

resetpoisonconc.o: resetpoisonconc.c header.h locations.h
	$(CC) $(CFLAGS) -c resetpoisonconc.c

resettimer.o: resettimer.c header.h
	$(CC) $(CFLAGS) -c resettimer.c

resizedynsrc.o: resizedynsrc.c header.h locations.h
	$(CC) $(CFLAGS) -c resizedynsrc.c

resizefissionsrc.o: resizefissionsrc.c header.h locations.h
	$(CC) $(CFLAGS) -c resizefissionsrc.c

responsefunction.o: responsefunction.c header.h locations.h
	$(CC) $(CFLAGS) -c responsefunction.c

riacycle.o: riacycle.c header.h locations.h
	$(CC) $(CFLAGS) -c riacycle.c

rroutput.o: rroutput.c header.h locations.h
	$(CC) $(CFLAGS) -c rroutput.c

runfinix.o: runfinix.c header.h locations.h
	$(CC) $(CFLAGS) -c runfinix.c

sabscattering.o: sabscattering.c header.h locations.h
	$(CC) $(CFLAGS) -c sabscattering.c

sampledelnu.o: sampledelnu.c header.h locations.h
	$(CC) $(CFLAGS) -c sampledelnu.c

sampleendflaw.o: sampleendflaw.c header.h locations.h
	$(CC) $(CFLAGS) -c sampleendflaw.c

sampleifcdata.o: sampleifcdata.c header.h locations.h
	$(CC) $(CFLAGS) -c sampleifcdata.c

samplemeshdelnu.o: samplemeshdelnu.c header.h locations.h
	$(CC) $(CFLAGS) -c samplemeshdelnu.c

samplemu.o: samplemu.c header.h locations.h
	$(CC) $(CFLAGS) -c samplemu.c

samplenu.o: samplenu.c header.h locations.h
	$(CC) $(CFLAGS) -c samplenu.c

sampleplasmasrc.o: sampleplasmasrc.c header.h locations.h
	$(CC) $(CFLAGS) -c sampleplasmasrc.c

samplepointdelnu.o: samplepointdelnu.c header.h locations.h
	$(CC) $(CFLAGS) -c samplepointdelnu.c

sampleprecursorgroup.o: sampleprecursorgroup.c header.h locations.h
	$(CC) $(CFLAGS) -c sampleprecursorgroup.c

sampleptable.o: sampleptable.c header.h locations.h
	$(CC) $(CFLAGS) -c sampleptable.c

samplereaction.o: samplereaction.c header.h locations.h
	$(CC) $(CFLAGS) -c samplereaction.c

samplesrcpoint.o: samplesrcpoint.c header.h locations.h
	$(CC) $(CFLAGS) -c samplesrcpoint.c

sampletabular.o: sampletabular.c header.h locations.h
	$(CC) $(CFLAGS) -c sampletabular.c

senscollision.o: senscollision.c header.h locations.h
	$(CC) $(CFLAGS) -c senscollision.c

scalarprod.o: scalarprod.c header.h locations.h
	$(CC) $(CFLAGS) -c scalarprod.c

schurfactorization.o: schurfactorization.c header.h locations.h
	$(CC) $(CFLAGS) -c schurfactorization.c

score.o: score.c header.h locations.h
	$(CC) $(CFLAGS) -c score.c

scorealb.o: scorealb.c header.h locations.h
	$(CC) $(CFLAGS) -c scorealb.c

scoreactidet.o: scoreactidet.c header.h locations.h
	$(CC) $(CFLAGS) -c scoreactidet.c

scoredf.o: scoredf.c header.h locations.h
	$(CC) $(CFLAGS) -c scoredf.c

scoreotfburn.o: scoreotfburn.c header.h locations.h
	$(CC) $(CFLAGS) -c scoreotfburn.c

scorepb.o: scorepb.c header.h locations.h
	$(CC) $(CFLAGS) -c scorepb.c

scoreperturbations.o: scoreperturbations.c header.h locations.h
	$(CC) $(CFLAGS) -c scoreperturbations.c

scorephotonheat.o: scorephotonheat.c header.h locations.h
	$(CC) $(CFLAGS) -c scorephotonheat.c

scoredirectsensitivity.o: scoredirectsensitivity.c header.h locations.h
	$(CC) $(CFLAGS) -c scoredirectsensitivity.c

scorepinpower.o: scorepinpower.c header.h locations.h
	$(CC) $(CFLAGS) -c scorepinpower.c

scorepoison.o: scorepoison.c header.h locations.h
	$(CC) $(CFLAGS) -c scorepoison.c

scorecapture.o: scorecapture.c header.h locations.h
	$(CC) $(CFLAGS) -c scorecapture.c

scorecpd.o: scorecpd.c header.h locations.h
	$(CC) $(CFLAGS) -c scorecpd.c

scorecmm.o: scorecmm.c header.h locations.h
	$(CC) $(CFLAGS) -c scorecmm.c

scorefet.o: scorefet.c header.h locations.h
	$(CC) $(CFLAGS) -c scorefet.c

scorefission.o: scorefission.c header.h locations.h
	$(CC) $(CFLAGS) -c scorefission.c

scoregc.o: scoregc.c header.h locations.h
	$(CC) $(CFLAGS) -c scoregc.c

scoreicmcol.o: scoreicmcol.c header.h locations.h
	$(CC) $(CFLAGS) -c scoreicmcol.c

scoreicmtrk.o: scoreicmtrk.c header.h locations.h
	$(CC) $(CFLAGS) -c scoreicmtrk.c

scoreinterfaceflux.o: scoreinterfaceflux.c header.h locations.h
	$(CC) $(CFLAGS) -c scoreinterfaceflux.c

scoreinterfacepower.o: scoreinterfacepower.c header.h locations.h
	$(CC) $(CFLAGS) -c scoreinterfacepower.c

scoreiternucs.o: scoreiternucs.c header.h locations.h
	$(CC) $(CFLAGS) -c scoreiternucs.c

scoremesh.o: scoremesh.c header.h locations.h
	$(CC) $(CFLAGS) -c scoremesh.c

scoremicrodep.o: scoremicrodep.c header.h locations.h
	$(CC) $(CFLAGS) -c scoremicrodep.c

scorermxcurr.o: scorermxcurr.c header.h locations.h
	$(CC) $(CFLAGS) -c scorermxcurr.c

scorermxmfp.o: scorermxmfp.c header.h locations.h
	$(CC) $(CFLAGS) -c scorermxmfp.c

scorermxresp.o: scorermxresp.c header.h locations.h
	$(CC) $(CFLAGS) -c scorermxresp.c

scorermxsrc.o: scorermxsrc.c header.h locations.h
	$(CC) $(CFLAGS) -c scorermxsrc.c

scorescattering.o: scorescattering.c header.h locations.h
	$(CC) $(CFLAGS) -c scorescattering.c

scorescaresp.o: scorescaresp.c header.h locations.h
	$(CC) $(CFLAGS) -c scorescaresp.c

scoresenslabel.o: scoresenslabel.c header.h locations.h
	$(CC) $(CFLAGS) -c scoresenslabel.c

scoresurf.o: scoresurf.c header.h locations.h
	$(CC) $(CFLAGS) -c scoresurf.c

scoretimeconstants.o: scoretimeconstants.c header.h locations.h
	$(CC) $(CFLAGS) -c scoretimeconstants.c

scoretimesource.o: scoretimesource.c header.h locations.h
	$(CC) $(CFLAGS) -c scoretimesource.c

scoretransmuxs.o: scoretransmuxs.c header.h locations.h
	$(CC) $(CFLAGS) -c scoretransmuxs.c

scoreufs.o: scoreufs.c header.h locations.h
	$(CC) $(CFLAGS) -c scoreufs.c

searcharray.o: searcharray.c header.h
	$(CC) $(CFLAGS) -c searcharray.c

seekarray.o: seekarray.c header.h
	$(CC) $(CFLAGS) -c seekarray.c

seeklist.o: seeklist.c header.h locations.h
	$(CC) $(CFLAGS) -c seeklist.c

seekliststr.o: seekliststr.c header.h locations.h
	$(CC) $(CFLAGS) -c seekliststr.c

sendddparticle.o: sendddparticle.c header.h locations.h
	$(CC) $(CFLAGS) -c sendddparticle.c

sendddparticlecounter.o: sendddparticlecounter.c header.h locations.h
	$(CC) $(CFLAGS) -c sendddparticlecounter.c

sendifcinputtemplates.o: sendifcinputtemplates.c header.h locations.h krakenmeshtypes.h
	$(CC) $(CFLAGS) -c sendifcinputtemplates.c

sendifcmesh.o: sendifcmesh.c header.h locations.h krakenmeshtypes.h
	$(CC) $(CFLAGS) -c sendifcmesh.c

sendifcoutputdata.o: sendifcoutputdata.c header.h locations.h
	$(CC) $(CFLAGS) -c sendifcoutputdata.c

sendifcoutputtemplates.o: sendifcoutputtemplates.c header.h locations.h krakenmeshtypes.h
	$(CC) $(CFLAGS) -c sendifcoutputtemplates.c

sensbufval.o: sensbufval.c header.h locations.h
	$(CC) $(CFLAGS) -c sensbufval.c

sensitivityoutput.o: sensitivityoutput.c header.h locations.h
	$(CC) $(CFLAGS) -c sensitivityoutput.c

separatebanuc.o: separatebanuc.c header.h locations.h
	$(CC) $(CFLAGS) -c separatebanuc.c

setcipopulation.o: setcipopulation.c header.h locations.h
	$(CC) $(CFLAGS) -c setcipopulation.c

setcoefcalc.o: setcoefcalc.c header.h locations.h
	$(CC) $(CFLAGS) -c setcoefcalc.c

setddidsimple.o: setddidsimple.c header.h locations.h
	$(CC) $(CFLAGS) -c setddidsimple.c

setdepstepsize.o: setdepstepsize.c header.h locations.h
	$(CC) $(CFLAGS) -c setdepstepsize.c

setdetflags.o: setdetflags.c header.h locations.h
	$(CC) $(CFLAGS) -c setdetflags.c

setdirectpointers.o: setdirectpointers.c header.h locations.h
	$(CC) $(CFLAGS) -c setdirectpointers.c

setfisse.o: setfisse.c header.h locations.h
	$(CC) $(CFLAGS) -c setfisse.c

sethisvcalc.o: sethisvcalc.c header.h locations.h
	$(CC) $(CFLAGS) -c sethisvcalc.c

setnormalization.o: setnormalization.c header.h locations.h
	$(CC) $(CFLAGS) -c setnormalization.c

setifctmslimits.o: setifctmslimits.c header.h locations.h
	$(CC) $(CFLAGS) -c setifctmslimits.c

setoption.o: setoption.c header.h locations.h
	$(CC) $(CFLAGS) -c setoption.c

setoptimization.o: setoptimization.c header.h locations.h
	$(CC) $(CFLAGS) -c setoptimization.c

setstlmeshpointers.o: setstlmeshpointers.c header.h locations.h
	$(CC) $(CFLAGS) -c setstlmeshpointers.c

setpathlevels.o: setpathlevels.c header.h locations.h
	$(CC) $(CFLAGS) -c setpathlevels.c

setprecursorgroups.o: setprecursorgroups.c header.h locations.h
	$(CC) $(CFLAGS) -c setprecursorgroups.c

shareinputdata.o: shareinputdata.c header.h locations.h
	$(CC) $(CFLAGS) -c shareinputdata.c

shuntingyard.o: shuntingyard.c header.h locations.h
	$(CC) $(CFLAGS) -c shuntingyard.c

signalexternal.o: signalexternal.c header.h locations.h
	$(CC) $(CFLAGS) -c signalexternal.c

signalhandler.o: signalhandler.c header.h locations.h
	$(CC) $(CFLAGS) -c signalhandler.c

socketreceivedouble.o: socketreceivedouble.c header.h locations.h
	$(CC) $(CFLAGS) -c socketreceivedouble.c

socketreceive.o: socketreceive.c header.h locations.h
	$(CC) $(CFLAGS) -c socketreceive.c

socketreceivelong.o: socketreceivelong.c header.h locations.h
	$(CC) $(CFLAGS) -c socketreceivelong.c

socketreceivestring.o: socketreceivestring.c header.h locations.h
	$(CC) $(CFLAGS) -c socketreceivestring.c

socketsenddouble.o: socketsenddouble.c header.h locations.h
	$(CC) $(CFLAGS) -c socketsenddouble.c

socketsend.o: socketsend.c header.h locations.h
	$(CC) $(CFLAGS) -c socketsend.c

socketsendlong.o: socketsendlong.c header.h locations.h
	$(CC) $(CFLAGS) -c socketsendlong.c

socketsendstring.o: socketsendstring.c header.h locations.h
	$(CC) $(CFLAGS) -c socketsendstring.c

solveotfburn.o: solveotfburn.c header.h locations.h
	$(CC) $(CFLAGS) -c solveotfburn.c

solvermx.o: solvermx.c header.h locations.h
	$(CC) $(CFLAGS) -c solvermx.c

solvermxcrit.o: solvermxcrit.c header.h locations.h
	$(CC) $(CFLAGS) -c solvermxcrit.c

sortall.o: sortall.c header.h locations.h
	$(CC) $(CFLAGS) -c sortall.c

sortarray.o: sortarray.c header.h
	$(CC) $(CFLAGS) -c sortarray.c

sortlist.o: sortlist.c header.h locations.h
	$(CC) $(CFLAGS) -c sortlist.c

speed.o: speed.c header.h locations.h
	$(CC) $(CFLAGS) -c speed.c

splitlist.o: splitlist.c header.h locations.h
	$(CC) $(CFLAGS) -c splitlist.c

srcdet.o: srcdet.c header.h locations.h
	$(CC) $(CFLAGS) -c srcdet.c

srcroulette.o: srcroulette.c header.h locations.h
	$(CC) $(CFLAGS) -c srcroulette.c

starttimer.o: starttimer.c header.h
	$(CC) $(CFLAGS) -c starttimer.c

statbin.o: statbin.c header.h locations.h
	$(CC) $(CFLAGS) -c statbin.c

statsum.o: statsum.c header.h locations.h
	$(CC) $(CFLAGS) -c statsum.c

stattests.o: stattests.c header.h locations.h
	$(CC) $(CFLAGS) -c stattests.c

stdcomp.o: stdcomp.c header.h natural_elements.h
	$(CC) $(CFLAGS) -c stdcomp.c

stddev.o: stddev.c header.h locations.h
	$(CC) $(CFLAGS) -c stddev.c

stlfacetdistance.o: stlfacetdistance.c header.h locations.h
	$(CC) $(CFLAGS) -c stlfacetdistance.c

stlmatfinder.o: stlmatfinder.c header.h locations.h
	$(CC) $(CFLAGS) -c stlmatfinder.c

stlraytest.o: stlraytest.c header.h locations.h
	$(CC) $(CFLAGS) -c stlraytest.c

stopatboundary.o: stopatboundary.c header.h locations.h
	$(CC) $(CFLAGS) -c stopatboundary.c

stopatwwbound.o: stopatwwbound.c header.h locations.h
	$(CC) $(CFLAGS) -c stopatwwbound.c

stopci.o: stopci.c header.h locations.h
	$(CC) $(CFLAGS) -c stopci.c

stopcciter.o: stopcciter.c header.h locations.h
	$(CC) $(CFLAGS) -c stopcciter.c

stopsurfdetflg.o: stopsurfdetflg.c header.h locations.h
	$(CC) $(CFLAGS) -c stopsurfdetflg.c

stoptimer.o: stoptimer.c header.h
	$(CC) $(CFLAGS) -c stoptimer.c

storecomposition.o: storecomposition.c header.h locations.h
	$(CC) $(CFLAGS) -c storecomposition.c

storehistorypoint.o: storehistorypoint.c header.h locations.h
	$(CC) $(CFLAGS) -c storehistorypoint.c

storermxevent.o: storermxevent.c header.h locations.h
	$(CC) $(CFLAGS) -c storermxevent.c

storesensevent.o: storesensevent.c header.h locations.h
	$(CC) $(CFLAGS) -c storesensevent.c

storetransmuxs.o: storetransmuxs.c header.h locations.h
	$(CC) $(CFLAGS) -c storetransmuxs.c

storevaluepair.o: storevaluepair.c header.h locations.h
	$(CC) $(CFLAGS) -c storevaluepair.c

sumdivcompositions.o: sumdivcompositions.c header.h locations.h
	$(CC) $(CFLAGS) -c sumdivcompositions.c

sumprivatedata.o: sumprivatedata.c header.h locations.h
	$(CC) $(CFLAGS) -c sumprivatedata.c

sumprivateres.o: sumprivateres.c header.h locations.h
	$(CC) $(CFLAGS) -c sumprivateres.c

sumtotxs.o: sumtotxs.c header.h locations.h
	$(CC) $(CFLAGS) -c sumtotxs.c

superdet.o: superdet.c header.h locations.h
	$(CC) $(CFLAGS) -c superdet.c

surfacedistance.o: surfacedistance.c header.h locations.h
	$(CC) $(CFLAGS) -c surfacedistance.c

supervisedcycle.o: supervisedcycle.c header.h locations.h krakensignals.h
	$(CC) $(CFLAGS) -c supervisedcycle.c

surfacenormal.o: surfacenormal.c header.h locations.h
	$(CC) $(CFLAGS) -c surfacenormal.c

surfacesrc.o: surfacesrc.c header.h locations.h
	$(CC) $(CFLAGS) -c surfacesrc.c

surfacevol.o: surfacevol.c header.h locations.h
	$(CC) $(CFLAGS) -c surfacevol.c

swapitems.o: swapitems.c header.h locations.h
	$(CC) $(CFLAGS) -c swapitems.c

swapuniverses.o: swapuniverses.c header.h locations.h
	$(CC) $(CFLAGS) -c swapuniverses.c

symboliclu.o: symboliclu.c header.h
	$(CC) $(CFLAGS) -c symboliclu.c

symmetryboundary.o: symmetryboundary.c header.h locations.h
	$(CC) $(CFLAGS) -c symmetryboundary.c

systemstat.o: systemstat.c header.h locations.h
	$(CC) $(CFLAGS) -c systemstat.c

targetvelocity.o: targetvelocity.c header.h locations.h
	$(CC) $(CFLAGS) -c targetvelocity.c

tmpmajorants.o: tmpmajorants.c header.h locations.h
	$(CC) $(CFLAGS) -c tmpmajorants.c

testasciifile.o: testasciifile.c header.h
	$(CC) $(CFLAGS) -c testasciifile.c

testdosfile.o: testdosfile.c header.h
	$(CC) $(CFLAGS) -c testdosfile.c

testfacetoverlap.o: testfacetoverlap.c header.h locations.h
	$(CC) $(CFLAGS) -c testfacetoverlap.c

testhisvbreak.o: testhisvbreak.c header.h locations.h
	$(CC) $(CFLAGS) -c testhisvbreak.c

testparam.o: testparam.c header.h
	$(CC) $(CFLAGS) -c testparam.c

testsurface.o: testsurface.c header.h locations.h
	$(CC) $(CFLAGS) -c testsurface.c

teststlgeometry.o: teststlgeometry.c header.h locations.h
	$(CC) $(CFLAGS) -c teststlgeometry.c

teststlsolids.o: teststlsolids.c header.h locations.h
	$(CC) $(CFLAGS) -c teststlsolids.c

testunisym.o: testunisym.c header.h locations.h
	$(CC) $(CFLAGS) -c testunisym.c

testvaluepair.o: testvaluepair.c header.h locations.h
	$(CC) $(CFLAGS) -c testvaluepair.c

testxs.o: testxs.c header.h locations.h
	$(CC) $(CFLAGS) -c testxs.c

tetputboundingbox.o: tetputboundingbox.c header.h locations.h
	$(CC) $(CFLAGS) -c tetputboundingbox.c

tetravol.o: tetravol.c header.h locations.h
	$(CC) $(CFLAGS) -c tetravol.c

thingrid.o: thingrid.c header.h
	$(CC) $(CFLAGS) -c thingrid.c

thompsonscattering.o: thompsonscattering.c header.h locations.h
	$(CC) $(CFLAGS) -c thompsonscattering.c

timecutoff.o: timecutoff.c header.h locations.h
	$(CC) $(CFLAGS) -c timecutoff.c

timeintervalstr.o: timeintervalstr.c header.h
	$(CC) $(CFLAGS) -c timeintervalstr.c

timerval.o: timerval.c header.h
	$(CC) $(CFLAGS) -c timerval.c

timercpuval.o: timercpuval.c header.h
	$(CC) $(CFLAGS) -c timercpuval.c

timestamp.o: timestamp.c header.h
	$(CC) $(CFLAGS) -c timestamp.c

timestr.o: timestr.c header.h
	$(CC) $(CFLAGS) -c timestr.c

tobank.o: tobank.c header.h locations.h
	$(CC) $(CFLAGS) -c tobank.c

tocommonque.o: tocommonque.c header.h locations.h
	$(CC) $(CFLAGS) -c tocommonque.c

tolimbo.o: tolimbo.c header.h locations.h
	$(CC) $(CFLAGS) -c tolimbo.c

toque.o: toque.c header.h locations.h
	$(CC) $(CFLAGS) -c toque.c

torusdis.o: torusdis.c header.h locations.h
	$(CC) $(CFLAGS) -c torusdis.c

tostack.o: tostack.c header.h locations.h
	$(CC) $(CFLAGS) -c tostack.c

tostore.o: tostore.c header.h locations.h
	$(CC) $(CFLAGS) -c tostore.c

totxs.o: totxs.c header.h locations.h
	$(CC) $(CFLAGS) -c totxs.c

trackfile.o: trackfile.c header.h locations.h
	$(CC) $(CFLAGS) -c trackfile.c

tracking.o: tracking.c header.h locations.h
	$(CC) $(CFLAGS) -c tracking.c

trackingerror.o: trackingerror.c header.h locations.h
	$(CC) $(CFLAGS) -c trackingerror.c

trackmode.o: trackmode.c header.h locations.h
	$(CC) $(CFLAGS) -c trackmode.c

transportcycle.o: transportcycle.c header.h locations.h
	$(CC) $(CFLAGS) -c transportcycle.c

transportcorrection.o: transportcorrection.c header.h locations.h
	$(CC) $(CFLAGS) -c transportcorrection.c

trapz.o: trapz.c header.h locations.h
	$(CC) $(CFLAGS) -c trapz.c

trapzreal.o: trapzreal.c header.h
	$(CC) $(CFLAGS) -c trapzreal.c

truncate.o: truncate.c header.h
	$(CC) $(CFLAGS) -c truncate.c

tta.o: tta.c header.h locations.h
	$(CC) $(CFLAGS) -c tta.c

ttachain.o: ttachain.c header.h
	$(CC) $(CFLAGS) -c ttachain.c

ttaloop.o: ttaloop.c header.h locations.h
	$(CC) $(CFLAGS) -c ttaloop.c

ttb.o: ttb.c header.h locations.h
	$(CC) $(CFLAGS) -c ttb.c

ufsfactor.o: ufsfactor.c header.h locations.h
	$(CC) $(CFLAGS) -c ufsfactor.c

unionizegrid.o: unionizegrid.c header.h locations.h
	$(CC) $(CFLAGS) -c unionizegrid.c

unisym.o: unisym.c header.h locations.h
	$(CC) $(CFLAGS) -c unisym.c

universeboundaries.o: universeboundaries.c header.h locations.h
	$(CC) $(CFLAGS) -c universeboundaries.c

updatecistop.o: updatecistop.c header.h locations.h
	$(CC) $(CFLAGS) -c updatecistop.c

updatefinixifc.o: updatefinixifc.c header.h locations.h
	$(CC) $(CFLAGS) -c updatefinixifc.c

updatefinixpower.o: updatefinixpower.c header.h locations.h
	$(CC) $(CFLAGS) -c updatefinixpower.c

updateifcdensmax.o: updateifcdensmax.c header.h locations.h
	$(CC) $(CFLAGS) -c updateifcdensmax.c

updateifctempminmax.o: updateifctempminmax.c header.h locations.h
	$(CC) $(CFLAGS) -c updateifctempminmax.c

updatemicrodens.o: updatemicrodens.c header.h locations.h
	$(CC) $(CFLAGS) -c updatemicrodens.c

updatermxwgt.o: updatermxwgt.c header.h locations.h
	$(CC) $(CFLAGS) -c updatermxwgt.c

uresdilumicroxs.o: uresdilumicroxs.c header.h locations.h
	$(CC) $(CFLAGS) -c uresdilumicroxs.c

uresfactor.o: uresfactor.c header.h locations.h
	$(CC) $(CFLAGS) -c uresfactor.c

userifc.o: userifc.c header.h locations.h
	$(CC) $(CFLAGS) -c userifc.c

useriter.o: useriter.c header.h locations.h
	$(CC) $(CFLAGS) -c useriter.c

usersrc.o: usersrc.c header.h locations.h
	$(CC) $(CFLAGS) -c usersrc.c

usersurf.o: usersurf.c header.h locations.h
	$(CC) $(CFLAGS) -c usersurf.c

usertransmucut.o: usertransmucut.c header.h locations.h
	$(CC) $(CFLAGS) -c usertransmucut.c

valuepairidx.o: valuepairidx.c header.h locations.h
	$(CC) $(CFLAGS) -c valuepairidx.c

valuepairval.o: valuepairval.c header.h locations.h
	$(CC) $(CFLAGS) -c valuepairval.c

vectornorm.o: vectornorm.c header.h locations.h
	$(CC) $(CFLAGS) -c vectornorm.c

vectorprod.o: vectorprod.c header.h locations.h
	$(CC) $(CFLAGS) -c vectorprod.c

virtgcucolflags.o: virtgcucolflags.c header.h locations.h
	$(CC) $(CFLAGS) -c virtgcucolflags.c

volumesmc.o: volumesmc.c header.h locations.h
	$(CC) $(CFLAGS) -c volumesmc.c

vrcycle.o: vrcycle.c header.h locations.h
	$(CC) $(CFLAGS) -c vrcycle.c

walkeralias.o: walkeralias.c header.h locations.h
	$(CC) $(CFLAGS) -c walkeralias.c

warn.o: warn.c header.h locations.h
	$(CC) $(CFLAGS) -c warn.c

weightwindow.o: weightwindow.c header.h locations.h
	$(CC) $(CFLAGS) -c weightwindow.c

whereami.o: whereami.c header.h locations.h
	$(CC) $(CFLAGS) -c whereami.c

workarray.o: workarray.c header.h locations.h
	$(CC) $(CFLAGS) -c workarray.c

writecimomfluxes.o: writecimomfluxes.c header.h locations.h
	$(CC) $(CFLAGS) -c writecimomfluxes.c

writedepfile.o: writedepfile.c header.h locations.h
	$(CC) $(CFLAGS) -c writedepfile.c

writedynsrc.o: writedynsrc.c header.h locations.h
	$(CC) $(CFLAGS) -c writedynsrc.c

writefinixifc.o: writefinixifc.c header.h locations.h
	$(CC) $(CFLAGS) -c writefinixifc.c

writefinixinputfile.o: writefinixinputfile.c header.h locations.h
	$(CC) $(CFLAGS) -c writefinixinputfile.c

writetetmeshtogeo.o: writetetmeshtogeo.c header.h locations.h
	$(CC) $(CFLAGS) -c writetetmeshtogeo.c

writeicmdata.o: writeicmdata.c header.h locations.h
	$(CC) $(CFLAGS) -c writeicmdata.c

writesensresult.o: writesensresult.c header.h locations.h
	$(CC) $(CFLAGS) -c writesensresult.c

writesourcefile.o: writesourcefile.c header.h locations.h
	$(CC) $(CFLAGS) -c writesourcefile.c

writestlmesh.o: writestlmesh.c header.h locations.h
	$(CC) $(CFLAGS) -c writestlmesh.c

writeumshtostl.o: writeumshtostl.c header.h locations.h
	$(CC) $(CFLAGS) -c writeumshtostl.c

writewwmesh.o: writewwmesh.c header.h locations.h
	$(CC) $(CFLAGS) -c writewwmesh.c

wwdis.o: wwdis.c header.h locations.h
	$(CC) $(CFLAGS) -c wwdis.c

wwimportance.o: wwimportance.c header.h locations.h
	$(CC) $(CFLAGS) -c wwimportance.c

wwinsrc.o: wwinsrc.c header.h locations.h
	$(CC) $(CFLAGS) -c wwinsrc.c

xsplotter.o: xsplotter.c header.h locations.h
	$(CC) $(CFLAGS) -c xsplotter.c

zaitoiso.o: zaitoiso.c header.h element_data.h
	$(CC) $(CFLAGS) -c zaitoiso.c

zdis.o: zdis.c header.h
	$(CC) $(CFLAGS) -c zdis.c

zonecount.o: zonecount.c header.h locations.h
	$(CC) $(CFLAGS) -c zonecount.c

###############################################################################

aux_functions.o: FINIX/aux_functions.c
	$(CC) $(CFLAGS) -c FINIX/aux_functions.c

clmech.o: FINIX/clmech.c
	$(CC) $(CFLAGS) -c FINIX/clmech.c

clmechaux.o: FINIX/clmechaux.c
	$(CC) $(CFLAGS) -c FINIX/clmechaux.c

clmechprop.o: FINIX/clmechprop.c
	$(CC) $(CFLAGS) -c FINIX/clmechprop.c

clthprop.o: FINIX/clthprop.c
	$(CC) $(CFLAGS) -c FINIX/clthprop.c

coolant.o: FINIX/coolant.c
	$(CC) $(CFLAGS) -c FINIX/coolant.c

finixdata.o: FINIX/finixdata.c
	$(CC) $(CFLAGS) -c FINIX/finixdata.c

finix_initialization.o: FINIX/finix_initialization.c
	$(CC) $(CFLAGS) -c FINIX/finix_initialization.c

finix_output.o: FINIX/finix_output.c
	$(CC) $(CFLAGS) -c FINIX/finix_output.c

fumech.o: FINIX/fumech.c
	$(CC) $(CFLAGS) -c FINIX/fumech.c

fumechprop.o: FINIX/fumechprop.c
	$(CC) $(CFLAGS) -c FINIX/fumechprop.c

futhprop.o: FINIX/futhprop.c
	$(CC) $(CFLAGS) -c FINIX/futhprop.c

gap.o: FINIX/gap.c
	$(CC) $(CFLAGS) -c FINIX/gap.c

heateq1d.o: FINIX/heateq1d.c
	$(CC) $(CFLAGS) -c FINIX/heateq1d.c

initial.o: FINIX/initial.c
	$(CC) $(CFLAGS) -c FINIX/initial.c

transient.o: FINIX/transient.c
	$(CC) $(CFLAGS) -c FINIX/transient.c

###############################################################################
