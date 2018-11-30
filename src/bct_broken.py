


# ***************************************************************************
#
#							   Bct:
#				de Bruijn graph based Short reads corrector for transcriptomic
#
#
#
# ***************************************************************************

# ############################################################################



debug_mode=0



def printCommand(cmd,pc=True):
	if (pc and debug_mode!=0):
		print(cmd,flush=True)

# get the platform
def getPlatform():
	if sys.platform == "linux" or sys.platform == "linux2":
		return "linux"
	elif sys.platform == "darwin":
		return "OSX"
	else:
		print("[ERROR] Bct is not compatible with Windows.")
		sys.exit(1);


# get the timestamp as string
def getTimestamp():
	return "[" + time.strftime("%H:%M:%S") + " " + time.strftime("%d/%m/%Y") + "] "



# check if reads files are present
def checkReadFiles(readfiles):
	if readfiles is None:
		return True
	allFilesAreOK = True
	#~ for file in readfiles:
	if not os.path.isfile(readfiles):
		print("[ERROR] File \""+file+"\" does not exist.")
		allFilesAreOK = False
	if not allFilesAreOK:
		dieToFatalError("One or more read files do not exist.")


# check if files written by Bct are present
def checkWrittenFiles(files):
	allFilesAreOK = True
	if not os.path.isfile(files):
		print("[ERROR] There was a problem writing \"" + files + "\".")
		allFilesAreOK = False
	if not allFilesAreOK:
		dieToFatalError("One or more files could not be written.")



# to return if an error makes the run impossible
def dieToFatalError (msg):
  print("[FATAL ERROR] " + msg)
  print("Try `Bct --help` for more information")
  sys.exit(1);


# launch subprocess
def subprocessLauncher(cmd, argstdout=None, argstderr=None,	 argstdin=None):
	args = shlex.split(cmd)
	p = subprocess.Popen(args, stdin = argstdin, stdout = argstdout, stderr = argstderr).communicate()
	return p

def printTime(msg, seconds):
	m, s = divmod(seconds, 60)
	h, m = divmod(m, 60)
	return msg + " %d:%02d:%02d" % (h, m, s)


def printWarningMsg(msg):
	print("[Warning] " + msg)


# ############################################################################
#			   graph generation with BCALM + BTRIM + BGREAT
# ############################################################################

def graphConstruction(Bct_MAIN, Bct_INSTDIR, OUT_DIR, fileBcalm, kmerSize, solidity, nb_cores, mappingEffort, missmatchAllowed,aSize,maximumOccurence,subsambleAnchor,alpha,low, OUT_LOG_FILES):
	try:
		inputBcalm=fileBcalm
		print("\n" + getTimestamp() + "--> Building the graph...",flush=True)
		os.chdir(OUT_LOG_FILES)
		logBcalm = "logBcalm"
		logBcalmToWrite = open(logBcalm, 'w')
		logTips = "logTips"
		logTipsToWrite = open(logTips, 'w')
		logBgreat = "logBgreat"
		logBgreatToWrite = open(logBgreat, 'w')
		#~ os.chdir(Bct_MAIN)
		os.chdir(OUT_DIR)
		coreUsed = "20" if nb_cores == 0 else str(nb_cores)

		if(os.path.isfile(OUT_DIR +"/dbg" + str(kmerSize)+".fa")):
			print("\t#Graph dbg" + str(kmerSize)+".fa: Already here ! Let us use it ", flush=True)
		else:
			print("\t#Graph  dbg" + str(kmerSize)+".fa: Construction... ", flush=True)
			# BCALM
			cmd=Bct_INSTDIR + "/bcalm -max-memory 10000 -in " + OUT_DIR + "/" + inputBcalm + " -kmer-size " + str(kmerSize) + " -abundance-min " + str(solidity) + " -out " + OUT_DIR + "/out " + " -nb-cores " + coreUsed

			printCommand( "\t\t"+cmd)
			p = subprocessLauncher(cmd, logBcalmToWrite, logBcalmToWrite)
			checkWrittenFiles(OUT_DIR + "/out.unitigs.fa")
			for filename in glob.glob(OUT_DIR + "/trashme*"):
				os.remove(filename)

			#  Graph Cleaning
			print("\t\t #Graph cleaning... ", flush=True)
			# BTRIM
			cmd=Bct_INSTDIR + "/btt -u out.unitigs.fa -k "+str(kmerSize)+" -t "+str(3*int(kmerSize-1))+" -T 5 -c "+coreUsed+" -o dbg"+str(kmerSize)+".fa -h  8 -a "+str(alpha)+" -l "+str(low)
			printCommand("\t\t\t"+cmd)
			p = subprocessLauncher(cmd, logTipsToWrite, logTipsToWrite)
			for filename in glob.glob(OUT_DIR + "/out.*"):
				os.remove(filename)
			os.remove("bankBcalm.txt")

		if(os.path.isfile(OUT_DIR +"/dbg" + str(kmerSize)+".fa")):
			# Read Mapping
			print("\t#Read mapping with BGREAT... ", flush=True)
			# BGREAT
			cmd=Bct_INSTDIR + "/bgreat -k " + str(kmerSize) + "  -u original_reads.fa -g dbg" + str(kmerSize) + ".fa -t " + coreUsed + " -a "+str(aSize)+" -o "+str(maximumOccurence)+" -i "+str(subsambleAnchor)+" -m "+str(missmatchAllowed)+" -c -O -f reads_corrected.fa -e "+str(mappingEffort)
			printCommand("\t\t"+cmd)
			p = subprocessLauncher(cmd, logBgreatToWrite, logBgreatToWrite)
			checkWrittenFiles(OUT_DIR + "/reads_corrected.fa")

		os.chdir(Bct_MAIN)

		print(getTimestamp() + "--> Done!")
		return {'kmerSize': kmerSize}
	except SystemExit:	# happens when checkWrittenFiles() returns an error
		sys.exit(1);
	except KeyboardInterrupt:
		sys.exit(1);
	except:
		print("Unexpected error during graph construction:", sys.exc_info()[0])
		dieToFatalError('')



# ############################################################################
#									Main
# ############################################################################
def main():

	wholeT = time.time()
	print("\n*** This is Bct - de Bruin graph based corrector for transcriptomic  ***\n")
	Bct_MAIN = os.path.dirname(os.path.realpath(__file__))
	print("Binaries are in: " + Bct_INSTDIR)

	# ========================================================================
	#						 Manage command line arguments
	# ========================================================================
	parser = argparse.ArgumentParser(description='Bct - De Bruijn graph based read corrector ',formatter_class=argparse.RawTextHelpFormatter)

	# ------------------------------------------------------------------------
	#							 Define allowed options
	# ------------------------------------------------------------------------
	parser.add_argument("-u", action="store", dest="single_readfiles",		type=str,					help="(MANDATORY) Input fasta read files. Several read files must be concatenated\n \n")
	parser.add_argument('-o', action="store", dest="out_dir",				type=str,	default=os.getcwd(),	help="Path to store the results (default = current directory)")
	parser.add_argument('-t', action="store", dest="nb_cores",				type=int,	default = 0,	help="Number of cores used (default max)")

	parser.add_argument('-k', action="store", dest="kSize",					type=int,	default = 31,	help="k-mer size (default 31)")
	parser.add_argument('-s', action="store", dest="min_cov",				type=int,	default = 2,	help="k-mer abundance threshold, k-mers present strictly less than this number of times in the dataset will be discarded (default 2)")
	parser.add_argument('-a', action="store", dest="relative_threshold",				type=int,	default = 10,	help="A path a time less covered than its alternative can be removed (default 10)")
	parser.add_argument('-l', action="store", dest="low_threshold",				type=int,	default = 5,	help="Suppicious patterns with a abundance inferior to l are removed (default 5)\n \n")
	#~ parser.add_argument('-S', action="store", dest="unitig_Coverage",				type=int,	default = 0,	help="unitig Coverage for  cleaning (default auto)\n")
	#~ parser.add_argument('-a', action="store", dest="aSize",	type=int,	default = 21,	help="an integer, Size of the anchor to use (default 21)")
	#~ parser.add_argument('-e', action="store", dest="mapping_Effort",				type=int,	default = 1000,	help="Anchors to test for mapping ")
	#~ parser.add_argument('-m', action="store", dest="missmatch_allowed",				type=int,	default = 10,	help="missmatch allowed in mapping (default 10)")
	parser.add_argument('-i', action="store", dest="subsamble_anchor",				type=int,	default = 1,	help="(ADVANCED) index one out of i anchors (default 1)")
	parser.add_argument('-n', action="store", dest="maximum_occurence",				type=int,	default = 8,	help="(ADVANCED) maximum occurence of an anchor (default 8)\n")
	parser.add_argument('-d', action="store", dest="DEBUG",				type=int,	default = 0,	help="(ADVANCED) Print command lines\n \n")


	#~ parser.add_argument('--version', action='version', version='%(prog)s 0.0.1')


	# ------------------------------------------------------------------------
	#				Parse and interpret command line arguments
	# ------------------------------------------------------------------------
	options = parser.parse_args()

	# ------------------------------------------------------------------------
	#				  Print command line
	# ------------------------------------------------------------------------
	print("The command line was: " + ' '.join(sys.argv))
	global debug_mode


	# ------------------------------------------------------------------------
	#				  Misc parameters
	# ------------------------------------------------------------------------
	kSize				= options.kSize
	min_cov				= options.min_cov
	aSize				= 21
	nb_cores			= options.nb_cores
	mappingEffort		= 1000
	#~ unitigCoverage		= options.unitig_Coverage
	missmatchAllowed		= 10
	alpha= options.relative_threshold
	low= options.low_threshold
	maximumOccurence		= options.maximum_occurence
	subsambleAnchor		= options.subsamble_anchor
	debug_mode		= options.DEBUG
	# ------------------------------------------------------------------------
	#				Create output dir and log files
	# ------------------------------------------------------------------------
	OUT_DIR = options.out_dir
	try:
		if not os.path.exists(OUT_DIR):
			os.mkdir(OUT_DIR)
		else:
			printWarningMsg(OUT_DIR + " directory already exists, Bct will use it.")

		outName = OUT_DIR.split("/")[-1]
		OUT_DIR = os.path.dirname(os.path.realpath(OUT_DIR)) + "/" + outName
		OUT_LOG_FILES = OUT_DIR + "/logs"
		if not os.path.exists(OUT_LOG_FILES):
			os.mkdir(OUT_LOG_FILES)
		parametersLog = open(OUT_DIR + "/ParametersUsed.txt", 'w');
		parametersLog.write("kSize:%s	k-mer_solidity:%s	aSize:%s	mapping_effort:%s	missmatch_allowed:%s maximum_occurence:%s subsample_anchor:%s\n " %(kSize, min_cov,aSize, mappingEffort,missmatchAllowed,maximumOccurence,subsambleAnchor))
		parametersLog.close()

		print("Results will be stored in: ", OUT_DIR)
	except:
		print("Could not write in out directory 1:", sys.exc_info()[0])
		dieToFatalError('')

	# ------------------------------------------------------------------------
	#				  Parse input read options
	# ------------------------------------------------------------------------
	try:
		bankBcalm = open(OUT_DIR + "/bankBcalm.txt", 'w');
	except:
		print("Could not write in out directory 2:", sys.exc_info()[0])

	# check if the given paired-end read files indeed exist
	paired_readfiles = None
	single_readfiles = None
	errorReadFile = 0
	paired_readfiles = None
	errorReadFile = 1

	# check if the given single-end read files indeed exist
	if options.single_readfiles:
		single_readfiles = ''.join(options.single_readfiles)
		try:
			single_readfiles = os.path.abspath(single_readfiles)
			checkReadFiles(options.single_readfiles)
			errorReadFile *= 0
		except:
			single_readfiles = None
			errorReadFile *= 1
	else:
		single_readfiles = None
		errorReadFile *= 1

	if errorReadFile:
		#~ parser.print_usage()
		parser.print_help()
		dieToFatalError("Bct requires at least a read file")

	bloocooArg = ""
	bgreatArg = ""
	paired = '' if paired_readfiles is None else str(paired_readfiles)
	single = '' if single_readfiles is None else str(single_readfiles)
	both = paired + "," + single



	if single_readfiles is not None and paired_readfiles is not None:  # paired end + single end
		fileCase = 3
		bankBcalm.write(OUT_DIR + "/reads_corrected1.fa\n" + OUT_DIR + "/reads_corrected2.fa\n")
	elif single_readfiles is None:	# paired end only
		fileCase = 1
		bankBcalm.write(OUT_DIR + "/original_reads.fa\n")
	else:  # single end only
		fileCase = 2
		bankBcalm.write(OUT_DIR + "/original_reads.fa\n")
	# bankBcalm.write(OUT_DIR + "lost_unitig.fa")
	bankBcalm.close()

	# ========================================================================
	#									RUN
	# ========================================================================


	# ------------------------------------------------------------------------
	#						   Kmer size selection
	# ------------------------------------------------------------------------
	t = time.time()
	os.chdir(OUT_DIR)
	cmd="ln -fs " + single_readfiles + " " + OUT_DIR + "/original_reads.fa"
	printCommand("\t\t\t"+cmd)
	p = subprocessLauncher(cmd)
	os.chdir(Bct_MAIN)



	# ------------------------------------------------------------------------
	#						   Graph construction and cleaning
	# ------------------------------------------------------------------------
	t = time.time()
	valuesGraph = graphConstruction(Bct_MAIN, Bct_INSTDIR, OUT_DIR, "bankBcalm.txt", kSize, min_cov, nb_cores, mappingEffort, missmatchAllowed,aSize,maximumOccurence,subsambleAnchor,alpha,low, OUT_LOG_FILES)
	print(printTime("Correction took: ", time.time() - t))






if __name__ == '__main__':
	main()
