
#include "miscRunner.hpp"
#include "elucidator/simulation.h"
#include <TwoBit.h>

#include "elucidator/seqToolsUtils.h"
// elucidator - A library for analyzing sequence data
// Copyright (C) 2012-2018 Nicholas Hathaway <nicholas.hathaway@umassmed.edu>,
//
// This file is part of elucidator.
//
// elucidator is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// elucidator is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with elucidator.  If not, see <http://www.gnu.org/licenses/>.
//

namespace njhseq {
static VecStr bioTerms = {"Absolute error","Adenine","Adenosine triphosphate","Adjacency edges","Adjacency list","Adjacent Nodes",
		"Affine gap penalty","Algorithm","Alignment","Alignment score","Alignment-based phylogeny","Allele","Allele frequency",
		"Alphabet","Alternative splicing","Alu repeat","Ambiguity symbol","Amino acid","Angstrom","Anticodon","Archaea","Array",
		"Array notation","Augmented string","Autosomal dominant disorder","Autosomal recessive disorder","Autosome","Average mass",
		"Average mass table","B Recognition Element","B-ion","Bacteriophage","Base pair","Basepair edges","Binary fission","Binary tree",
		"Binomial random variable","BioPython","Blending inheritance","BLOSUM62","Bonding Graph","Bowtie","Brute force","C-terminus",
		"Catalan numbers","Cell","Centimorgan","Central dogma of molecular biology","Character","Character table","Characterizable strings",
		"Child","Chimeric protein","Chromatin","Chromosome","Circular string","Clustal","Coding strand","Codon","Combination",
		"Combinatorics","Common logarithm","Common subsequence","Common substring","Common supersequence","Comparative genomics",
		"Complementary base","Complementary DNA","Complementary event","Complete graph","Complete spectrum","Condensation reaction",
		"Conditional probability","Connected graph","Consensus string","Consistent character table","Consistent characters",
		"Consistent matrix","Consistent quartet","Constant gap penalty","Contig","Convergence of characters","CpG island","CpG site",
		"CRISPR","Cut","Cycle","Cystic fibrosis","Cytosine","Dalton","Data structure","de Bruijn graph","Decreasing permutation subsequence",
		"Degree","Dehydration reaction","Deoxyribose","Dependent random variables","Difference multiset","Dimer","Diploid cell",
		"Directed acyclic graph","Directed cycle","Directed Edge","Directed graph","Directed loop","Directed path",
		"Disjoint sets","Distance","Distance matrix","Distance-based phylogeny","Distinct","DNA","DNA codon table","DNA ligase",
		"DNA methylation","DNA polymerase","DNA replication","DNA string","DNAfull","Domain","Dominant allele","Double helix",
		"Dynamic programming","Edge","Edge weight","Edit alignment score","Edit distance","Edit operation","Element",
		"EMBOSS","Empty set","Endosymbiont","Entrez","Enumerate","Enzyme","Epigenetics",
		"Euclidean distance","Eukaryote","Eulerian cycle","European Bioinformatics Institute","European Molecular Biology Laboratory",
		"Exon","Exon skipping","Expected value","Factor","Failure array","FASTA format","FASTA package","FastQC","FASTX",
		"Fibonacci sequence","Fitting alignment","Founder effect","Frameshift mutation","Galaxy","Gamete","Gap","Gap deletion",
		"Gap insertion","Gap penalty","Gap symbol","GC-content","GC-poor","GC-rich","GCG","Gel electrophoresis","GenBank","Gene",
		"Gene coding region","Gene expression","Gene flow","Gene regulation","Gene symbol","Genetic carrier","Genetic code",
		"Genetic code variant","Genetic drift","Genetic equilibrium","Genetic fingerprinting","Genetic linkage",
		"Genetic recombination","Genetic string","Genome","Genome 10K Project","Genome assembly","Genome browser",
		"Genome rearrangement","Genome sequencing","Genomics","Genotype","Genotyping","Genus","GLAM2","Graph",
		"Graph algorithm","Graph theory","Graphviz","Guanine","Hairpin loop","Hamming distance","Haploid cell","Hardy-Weinberg principle","Head","Heredity","Heterozygous","Heuristic","Histone protein","Homodimer","Homologous","Homologous chromosomes","Homozygous","Homozygous dominant","Homozygous recessive","Hydration reaction","IDLE","Incident","Increasing permutation subsequence","Indel","Independent events","Independent random variables","Indicator random variable","Intensity","Internal Node","Intersection","Interwoven Strings","Intron","Inversion","Ion","Isotope","IUPAC notation","Junk DNA","k-fold substring","k-mer","k-mer composition","kbp","Knuth-Morris-Pratt algorithm","Leaf","Leslie matrix","Lexicographic order","Linear gap penalty","Linguistic complexity","Local alignment","Local alignment problem","Location","Longest common subsequence","Longest common substring","Loop","Macromolecule","Marker","Markov chain","Mass spectrometry","Mass spectrum","Matching","Matrix","Matrix multiplication","Maximal repeat","Maximum matching","Meiosis","MEME","Memorylessness","MeSH","Messenger RNA","Microsatellite","Minkowski Difference","Minkowski sum","Mismatch score","Mitochondrion","Mitosis","Modular arithmetic","Monoisotopic mass","Monoisotopic mass table","Monomer","Motif","Motif regular expression","Motzkin numbers","Multiple alignment","Multiple alignment score","Multiplicity","Multiset","Mutation","N statistic","N-glycosylation motif","N-terminus","National Center for Biotechnology Information","Natural selection","Neighbor","Newick format","NEXUS","Node","Noncrossing matching","Nontrivial character","Nontrivial split","Normalized discrepancy","Nucleic Acid","Nucleic acid primary structure","Nucleic acid secondary structure","Nucleic acid tertiary structure","Nucleobase","Nucleosome","Nucleotide","Nucleus","Numerical analysis","Oligonucleotide","Open reading frame","Optimal alignment","Optimize","Organelle","Outcome","Overlap alignment","Overlap graph","p-distance","PAM250","Parent","Parent mass","Parsimony","Partial character","Partial character table","Partial permutation","Partial split","Path","Path length","Pattern","Pattern matching","PDB","Peptide","Peptide bond","Perfect coverage","Perfect matching","Permutation","Permutation subsequence","Pfam","Phenotype","Phosphate","Phosphodiester bond","PHYLIP format","Phylogeny","PIR","Plasmid","Point mutation","Polymer","Polymerase chain reaction","Polymorphic trait","Polypeptide","Polyploid cell","Population bottleneck","Position","Precursor mRNA","Prefix","Prefix spectrum","Primer","Probabilistic event","Probability","Probability tree diagram","Profile matrix","Prokaryote","Promoter","PROSITE","Protein","Protein domain","Protein family","Protein isoform","Protein primary structure","Protein quaternary structure","Protein secondary structure","Protein string","Protein tertiary structure","Proteomics","Pseudoknot","PubMed","Punnett square","Purine","Pyrimidine","Python","Quartet","Quartet distance","Quartet-based phylogeny","Random string","Random variable","Rare-cutter enzyme","Read","Read coverage","Read generation","Reading frame","Recessive allele","Recognition sequence","Recurrence relation","RefSeq","Repeat","Repeated substring","Residue","Restriction digest","Restriction enzyme","Restriction map","Restriction site","Reversal","Reversal distance","Reversal sorting","Reverse complement","Reverse palindrome","Reverse transcription","Reversing substitution","Revert","Ribose","Ribosome","RNA","RNA codon table","RNA folding","RNA interference","RNA Polymerase","RNA splicing","RNA string","RNA transcription","Root","Rooted binary tree","Rooted Newick format","Rooted tree","SAM","ScanProsite","Scoring matrix","Semiglobal alignment","Sequence","Sequence logo","Set","Set complement","Set difference","Set notation","Set theory","Sex chromosome","Sex linkage","Shared peaks count","Shortest common supersequence","Sickle-cell anemia","Side chain","Signed permutation","Silent substitution","Simple cycle","Simplified spectrum","Single gene disorder","Single-nucleotide polymorphism","SMS 2","Somatic cell","Sorting","Sorting reversal","Spacer","Spectral convolution","Spectrum graph","Spliceosome","Split","Split distance","Start codon","Stop codon","Strand","String","String algorithm","String length","String weight","Submatrix","Suboptimal alignment","Subsequence","Subsequence indices","Subset","Substitution model","Substring","Subtree","Suffix","Suffix array","Suffix tree","Sugar","Sugar-phosphate backbone","Supersequence","Superstring","Symbol","Symbol weight","Synteny block","t-prefix","t-suffix","Tail","Tandem mass spectrometry","TATA Box","Taxon","Template strand","Text","Thymine","Trait","Transfer RNA","Transition","Transition matrix","Transition/transversion ratio","Translation","Transposon","Transversion","Tree","Tree of Life","Trie","Trimmomatic","Trivial character","Trivial split","Undirected graph","Uniform random variable","Union","UniProt","Unlabeled tree","Unrooted binary tree","Unrooted Newick format","Uracil","Valid basepair matching","Walk","Weighted alphabet","Weighted graph","Weighted string","Wild type","Wobble base pair","Wright-Fisher model","X Chromosome"};


miscRunner::miscRunner()
    : njh::progutils::ProgramRunner({addFunc("listAllFiles", listAllFiles, false),

                     addFunc("hangman", hangman, true),
                     addFunc("aminoAcidQuiz", aminoAcidQuiz, false),
                     addFunc("getSharedLines", getSharedLines, false),
										 addFunc("genInputForEstimateS", genInputForEstimateS, false),
										 addFunc("getHighestHapFrac", getHighestHapFrac, false),



										 addFunc("parseBamForForwardPerfectHits", parseBamForForwardPerfectHits, false),

										 addFunc("expandTableBySeparatingColumn", expandTableBySeparatingColumn, false),

										 addFunc("codeComparison", codeComparison, false),

										 addFunc("getSlidingQualityWindowMeans", getSlidingQualityWindowMeans, false),
										 addFunc("isFileEmpty", isFileEmpty, false),

										 addFunc("createConnectedHaplotypeNetwork", createConnectedHaplotypeNetwork, false),
										 addFunc("getLinkedInfoFromAdjList", getLinkedInfoFromAdjList, false),
										 addFunc("printSeqsColored", printSeqsColored, false),
										 addFunc("renameDownloadedGenome", renameDownloadedGenome, false),


										 addFunc("createSharedPathwaysFromReads", createSharedPathwaysFromReads, false),
										 addFunc("createSharedPathwaysFromContigs", createSharedPathwaysFromContigs, false),
										 addFunc("createSharedPathwaysFromRefSeqs", createSharedPathwaysFromRefSeqs, false),
										 addFunc("createSharedSubSegmentsFromRefSeqs", createSharedSubSegmentsFromRefSeqs, false),
										 },//
                    "misc") {}


//,,


int miscRunner::printSeqsColored(const njh::progutils::CmdArgs & inputCommands){
	OutOptions outOpts(bfs::path(""), ".fasta");
	seqSetUp setUp(inputCommands);
	setUp.processReadInNames();
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);
	OutputStream out(outOpts);
	SeqInput seqReader(setUp.pars_.ioOptions_);
	seqReader.openIn();
	seqInfo seq;
	while(seqReader.readNextRead(seq)){
		seq.outPutSeqAnsi(out);
	}
	return 0;
}


int miscRunner::isFileEmpty(const njh::progutils::CmdArgs & inputCommands){
	bfs::path fnp;
	seqSetUp setUp(inputCommands);
	setUp.setOption(fnp, "--fnp", "File path to file, can be gzipped and the contents of gzip file would be checked", true);
	setUp.finishSetUp(std::cout);

	std::stringstream ss;
	{
		InputStream in{InOptions(fnp)};
		//ss << in.rdbuf();
		std::string line;
		njh::files::crossPlatGetline(in,line);
		ss << line;
	}
	std::cout << njh::colorBool("" ==  ss.str()) << std::endl;
	return 0;
}

int miscRunner::parseBamForForwardPerfectHits(const njh::progutils::CmdArgs & inputCommands){
	uint32_t allowableErrors = 0;
	OutOptions outOpts(bfs::path(""));
	outOpts.outExtention_ = ".bed";
	bfs::path twoBitFnp = "";
	seqSetUp setUp(inputCommands);
	setUp.processReadInNames( { "--bam" }, true);
	setUp.setOption(allowableErrors, "--allowableErrors",
			"allowable errors in a motif element");
	setUp.setOption(twoBitFnp, "--twoBitFnp", "Two Bit Fnp", true);
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	comparison allowableErrorsComp;
	allowableErrorsComp.hqMismatches_ = allowableErrors;
	std::vector<std::shared_ptr<AlignmentResults>> alnResults = gatherMapResults(
			setUp.pars_.ioOptions_.firstName_, twoBitFnp, allowableErrorsComp);
	OutputStream out(outOpts);
	for(const auto & result : alnResults){
		if(!result->gRegion_.reverseSrand_){
			auto outRegion = result->gRegion_.genBedRecordCore();
			outRegion.name_ = result->bAln_.Name;
			out << outRegion.toDelimStr() << std::endl;
		}
	}

	return 0;
}


int miscRunner::genInputForEstimateS(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
	std::string filename = "";

	setUp.setOption(filename, "--file", "Filename of the output of SeekDeep processClusters", true);
	setUp.pars_.ioOptions_.firstName_ = filename;
	std::string cutOffs = "0.5,1,2.5,5,10";
	setUp.setOption(cutOffs, "--cutOffs", "A comma delimited string indicated the cut offs in percents");
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);

	table inTab(filename, "\t", true);
	std::vector<double> perCutOffs = vecStrToVecNum<double>(tokenizeString(cutOffs, ","));
	for(const auto & cut : perCutOffs){
		auto currentTab = inTab.extractNumColGreater("c_AveragedFrac", cut/100.0);
		std::set<std::string> popUidsSet;
		auto popStrs = currentTab.getColumn("h_popUID");
		for(const auto & str : popStrs){
			popUidsSet.emplace(str);
		}
		std::vector<std::string> popUids(popUidsSet.begin(), popUidsSet.end());

		auto tabSplit = currentTab.splitTableOnColumn("s_Sample");
		std::ofstream outFile;
		openTextFile(outFile, setUp.pars_.directoryName_ + estd::to_string(cut), ".txt", false, true);
		outFile << cut << std::endl;
		outFile << popUids.size() << " " << tabSplit.size() << std::endl;
		outFile << "\t" << vectorToString(popUids, "\t") << std::endl;;
		std::unordered_map<std::string,uint32_t> counts;
		for(const auto & split : tabSplit){
			outFile << split.first;
			auto currentPopUIDs = split.second.getColumn("h_popUID");
			for(const auto & pop : popUids){
				outFile << "\t";
				if(njh::in(pop, currentPopUIDs)){
					outFile << 1;
					++counts[pop];
				}else{
					outFile << 0;
				}
			}
			outFile << std::endl;
		}

		for(const auto & pop : popUids){
			outFile << "\t" << counts[pop];
		}
		outFile << std::endl;
	}
	return 0;
}

int miscRunner::getHighestHapFrac(const njh::progutils::CmdArgs & inputCommands){
	seqSetUp setUp(inputCommands);
	std::string filename = "";

	setUp.setOption(filename, "--file", "Filename of the output of SeekDeep processClusters", true);
	setUp.pars_.ioOptions_.firstName_ = filename;
	std::string getCols = "s_Sample,h_popUID,c_AveragedFrac";
	setUp.setOption(getCols, "--getCols", "A comma delimited string indicated the columns to extract");
	setUp.processDirectoryOutputName(true);
	setUp.finishSetUp(std::cout);

	table inTab(filename, "\t", true);
	VecStr cols = tokenizeString(getCols, ",");
	auto tabSplit = inTab.splitTableOnColumn("s_Sample");
	std::ofstream outFile;
	openTextFile(outFile, setUp.pars_.directoryName_, ".tab.txt", false, true);
	table finalTab(cols);
	for(const auto & split : tabSplit){
		auto cTab = split.second.extractByComp("c_clusterID", [](const std::string & str){
			return "0" == str;
		});
		finalTab.rbind(cTab.getColumns(cols), false);
	}
	finalTab.outPutContents(
			TableIOOpts(
					OutOptions(setUp.pars_.directoryName_ + "out", ".tab.txt", "tab",
							false, false, false), "\t", finalTab.hasHeader_));
	return 0;
}


int miscRunner::hangman(const njh::progutils::CmdArgs & inputCommands){
	std::string start = "     ___________\n     |/      |\n     |      \n     |      \n     |       \n     |       \n     |\n    _|___\n";
	std::string hang1 =		"     ___________\n     |/      |\n     |      (_)\n     |      \n     |       \n     |       \n     |\n    _|___\n";
	//std::string hang2 =		"     ___________\n     |/      |\n     |      (_)\n     |       |\n     |       \n     |       \n     |\n    _|___\n";
	std::string hang3 =		"     ___________\n     |/      |\n     |      (_)\n     |       |\n     |       |\n     |       \n     |\n    _|___\n";
	std::string hang4 =		"     ___________\n     |/      |\n     |      (_)\n     |      \\|\n     |       |\n     |       \n     |\n    _|___\n";
	std::string hang5 =		"     ___________\n     |/      |\n     |      (_)\n     |      \\|/\n     |       |\n     |       \n     |\n    _|___\n";
	std::string hang6 =		"     ___________\n     |/      |\n     |      (_)\n     |      \\|/\n     |       |\n     |      / \n     |\n    _|___\n";
	std::string hang7 =		"     ___________\n     |/      |\n     |      (_)\n     |      \\|/\n     |       |\n     |      / \\ \n     |\n    _|___\n";
	std::vector<std::string> hangMans ={start, hang1, hang3, hang4, hang5, hang6, hang7};
	uint32_t hangs = 0;
	std::vector<char> correctAttempts;
	std::vector<char> incorrectAttempts;
	//std::string ascaiiDna ="(==(     )==)\n `-.`. ,',-'\n    _,-'\"\n ,-',' `.`-.\n(==(     )==)\n `-.`. ,',-'\n    _,-'\"\n ,-',' `.`-.\n(==(     )==)\n `-.`. ,',-'\n    _,-'\"\n ,-',' `.`-.\n(==(     )==)\n `-.`. ,',-'\n    _,-'\"\n ,-',' `.`-.\n(==(     )==)";
	//std::cout << ascaiiDna << std::endl;
	randomGenerator gen;
	//std::string ans = "b";
	//std::string ans = "bioinformatics";
	std::string ans = bioTerms[gen.unifRand<uint32_t>(0, len(bioTerms))];

	std::string ansLower = stringToLowerReturn(ans);
	std::string currentAns = std::string(len(ans), '*');
	auto spaces = findOccurences(ans, " ");
	for(const auto & p : spaces){
		currentAns[p] = ans[p];
	}
	auto dashes = findOccurences(ans, "-");
	for(const auto & p : dashes){
		currentAns[p] = ans[p];
	}//"\033[" + atr + "m"
	std::string standard = "\033[38;5;231;48;5;16m";
njh::bashCT::resetAdd("");
	auto addFlashingColor = [](uint32_t colorCode){
		return njh::bashCT::flashing + njh::bashCT::addColor(colorCode);
	};
	//std::string win = boldText("                                  .''.\n       .''.      .        *''*    :_\\/_:     .\n      :_\\/_:   _\\(/_  .:.*_\\/_*   : /\\ :  .'.:.'.\n  .''.: /\\ :    /)\\   ':'* /\\ *  : '..'.  -=:o:=-\n :_\\/_:'.:::.  | ' *''*    * '.\\'/.'_\\(/_'.':'.'\n : /\\ : :::::  =  *_\\/_*     -= o =- /)\\    '  *\n  '..'  ':::' === * /\\ *     .'/.\\'.  ' ._____\n      *        |   *..*         :       |.   |' .---\"|\n        *      |     _           .--'|  ||   | _|    |\n        *      |  .-'|       __  |   |  |    ||      |\n     .-----.   |  |' |  ||  |  | |   |  |    ||      |\n ___'       ' /\"\\ |  '-.\"\".    '-'   '-.'    '`      |____\n   " + boldText("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~", "34") + "\n                      ~-~-~-~-~-~-~-~-~-~   /|\n          )      ~-~-~-~-~-~-~-~  /|~       /_|\\\n        _-H-__  -~-~-~-~-~-~     /_|\\    -~======-~\n~-\\XXXXXXXXXX/~     ~-~-~-~     /__|_\\ ~-~-~-~\n~-~-~-~-~-~    ~-~~-~-~-~-~    ========  ~-~-~-~\n", "45");
	//std::string win = boldText("                                  .''.\n       .''.      .        *''*    :_\\/_:     .\n      :_\\/_:   _\\(/_  .:.*_\\/_*   : /\\ :  .'.:.'.\n  .''.: /\\ :    /)\\   ':'* /\\ *  : '..'.  -=:o:=-\n :_\\/_:'.:::.  | ' *''*    * '.\\'/.'_\\(/_'.':'.'\n : /\\ : :::::  =  *_\\/_*     -= o =- /)\\    '  *\n  '..'  ':::' === * /\\ *     .'/.\\'.  ' ._____\n      *        |   *..*         :       |.   |' .---\"|\n        *      |     _           .--'|  ||   | _|    |\n        *      |  .-'|       __  |   |  |    ||      |\n     .-----.   |  |' |  ||  |  | |   |  |    ||      |\n ___'       ' /\"\\ |  '-.\"\".    '-'   '-.'    '`      |____\n   ", "45") + boldText(boldText("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~", "34"), "45") + boldText("\n                      " + boldText("~-~-~-~-~-~-~-~-~-~","34") + "   /|\n          )      " + boldText("~-~-~-~-~-~-~-~","34") + "  /|~       /_|\\\n        _-H-__  -~-~-~-~-~-~     /_|\\    -~======-~\n~-\\XXXXXXXXXX/~     ~-~-~-~     /__|_\\ ~-~-~-~\n~-~-~-~-~-~    ~-~~-~-~-~-~    ========  ~-~-~-~\n", "45");
	std::string line0  = standard;
	std::string line1  = "                                  " + addFlashingColor(197) + " .''." + njh::bashCT::resetAdd(standard);
	std::string line2  = "      " + addFlashingColor(202) + " .''.     " + addFlashingColor(83) + " .        " + addFlashingColor(128) + "*''*   " + addFlashingColor(197) + " :_\\/_:" + njh::bashCT::resetAdd(standard) + "     " + addFlashingColor(33) + "." + njh::bashCT::resetAdd(standard);
	std::string line3  = "      " + addFlashingColor(202) + ":_\\/_:   " + addFlashingColor(83) + "_\\(/_ " + addFlashingColor(51) + " .:." + addFlashingColor(128) + "*_\\/_*  " + addFlashingColor(197) + " : /\\ :" + njh::bashCT::resetAdd(standard) + "  " + addFlashingColor(33) + ".'.:.'." + njh::bashCT::resetAdd(standard);
	std::string line4  = "  " + addFlashingColor(46) + ".''." + addFlashingColor(202) + ": /\\ :    " + addFlashingColor(83) + "/)\\  " + addFlashingColor(51) + " ':'" + addFlashingColor(128) + "* /\\ *" + addFlashingColor(208) + " : " + addFlashingColor(197) + "'.''." + addFlashingColor(33) + "  -=:o:=-" + njh::bashCT::resetAdd(standard);
	std::string line5  = " " + addFlashingColor(46) + ":_\\/_:" + addFlashingColor(202) + "'" + addFlashingColor(51) + ".:::.  " + njh::bashCT::resetAdd(standard) + "| " + addFlashingColor(83) + "'" + addFlashingColor(201) + " *''*   " + addFlashingColor(128) + " *" + addFlashingColor(208) + " '.\\'/." + addFlashingColor(47) + "'_\\(/_'" + addFlashingColor(33) + ".':'.'" + njh::bashCT::resetAdd(standard);
	std::string line6  = " " + addFlashingColor(46) + ": /\\ :" + addFlashingColor(51) + " ::::: " + njh::bashCT::resetAdd(standard) + " =  " + addFlashingColor(201) + "*_\\/_*    " + addFlashingColor(208) + " -= o =-" + addFlashingColor(47) + " /)\\  " + addFlashingColor(33) + "  '  *" + njh::bashCT::resetAdd(standard);
	std::string line7  = "  " + addFlashingColor(46) + "'..' " + addFlashingColor(51) + " ':::'" + njh::bashCT::resetAdd(standard) + " === " + addFlashingColor(201) + "* /\\ *     " + addFlashingColor(208) + ".'/.\\'.  " + addFlashingColor(47) + "' " + njh::bashCT::resetAdd(standard) + "._____";
	std::string line8  = "   " + addFlashingColor(46) + "   * " + njh::bashCT::resetAdd(standard) + "       |   " + addFlashingColor(201) + "*..*         " + addFlashingColor(208) + ":      " + njh::bashCT::resetAdd(standard) + " |.   |' .---\"|";
	std::string line9  = "   " + addFlashingColor(46) + "     *  " + njh::bashCT::resetAdd(standard) + "    |     _           .--'|  ||   | _|    |";
	std::string line10 = "   " + addFlashingColor(46) + "     *  " + njh::bashCT::resetAdd(standard) + "    |  .-'|       __  |   |  |    ||      |";
	std::string line11 = "     .-----.   |  |' |  ||  |  | |   |  |    ||      |";
	std::string line12 = " ___'       ' /\"\\ |  '-.\"\".    '-'   '-.'    '`      |____";
	std::string line13 = "     " + njh::bashCT::addColor(21) +"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + njh::bashCT::resetAdd(standard);
	std::string line14 = "                      " + njh::bashCT::addColor(21) +"~-~-~-~-~-~-~-~-~-~" + njh::bashCT::resetAdd(standard) + "   /|";
	std::string line15 = "          )      " + njh::bashCT::addColor(21) +"~-~-~-~-~-~-~-~" + njh::bashCT::resetAdd(standard) + "  /|~       /_|\\";
	std::string line16 = "        _-H-__  " + njh::bashCT::addColor(21) +"-~-~-~-~-~-~" + njh::bashCT::resetAdd(standard) + "     /_|\\    -~======-~";
	std::string line17 = "~-\\XXXXXXXXXX/~    " + njh::bashCT::addColor(21) +" ~-~-~-~" + njh::bashCT::resetAdd(standard) + "     /__|_\\" + njh::bashCT::addColor(21) +" ~-~-~-~" + njh::bashCT::resetAdd(standard);
	std::string line18 = "" + njh::bashCT::addColor(21) +"~-~-~-~-~-~   " + njh::bashCT::addColor(21) +" " + njh::bashCT::addColor(21) +"~-~~-~-~-~-~" + njh::bashCT::resetAdd(standard) + "    ========  " + njh::bashCT::addColor(21) +"~-~-~-~" + njh::bashCT::reset;
	VecStr winLines = {line0,line1,line2,line3,line4,line5,line6,line7,line8,line9,line10,line11,
	line12,line13, line14, line15, line16, line17, line18};
	while(hangs < (len(hangMans) - 1) && currentAns != ans ){
		if(len(incorrectAttempts) + len(correctAttempts) > 0){
			for(uint32_t i = 0; i < 19; ++i){
				std::cout << "\r";
			}
			std::flush(std::cout);
		}
		std::cout << "current word" << std::endl;
		std::cout << currentAns << std::endl;
		std::cout << "attempts left" << std::endl;
		std::cout << len(hangMans) - 1 - hangs << std::endl;

		std::cout << njh::bashCT::bold
				<< njh::bashCT::green << "Current correct letters"
				<< njh::bashCT::reset << std::endl;
		njh::sort(correctAttempts);
		printVector(correctAttempts, ",");
		std::cout <<  njh::bashCT::bold
				<< njh::bashCT::red
				<< "Current incorrect letters"
				<< njh::bashCT::reset << std::endl;
		njh::sort(incorrectAttempts);
		printVector(incorrectAttempts, ",");
		std::string currentGuessStr = "";
		std::cout << hangMans[hangs] << std::endl;
		std::cout << "Enter guess" << std::endl;
		std::getline(std::cin, currentGuessStr);
		if(len(currentGuessStr) > 0){
			char currentGuess = currentGuessStr[0];
			currentGuess = tolower(currentGuess);
			if(vectorContains(incorrectAttempts, currentGuess) || vectorContains(correctAttempts, currentGuess)){
				std::cout << "Already guessed " << currentGuessStr[0] << std::endl;
			}else{
				if(ansLower.find(currentGuess) == std::string::npos){
					incorrectAttempts.emplace_back(currentGuess);
					++hangs;
				}else{
					correctAttempts.emplace_back(currentGuess);
					auto pos = findOccurences(ansLower, std::string(1,currentGuess));
					for(const auto & p : pos){
						currentAns[p] = ans[p];
					}
				}
			}

			//std::cout << hangMans[hangs] << std::endl;
			//std::cout << currentAns << std::endl;
		}
	}
	if(hangs < len(hangMans) - 1){
		std::cout << njh::bashCT::bold
				<< njh::bashCT::green
			<< "Winner!" << njh::bashCT::reset << std::endl;
		std::cout << njh::bashCT::bold
				<< ans << njh::bashCT::reset << std::endl;
		//std::cout << boldText(win, "40") << std::endl;
		printVector(winLines, "\n");
	}else{
		std::cout << njh::bashCT::bold
				<< njh::bashCT::red << "Lost!" << njh::bashCT::reset  << std::endl;
		std::cout << hangMans.back() << std::endl;
		std::cout << "Word was " << ans << std::endl;
	}
	return 0;
}









int miscRunner::listAllFiles(const njh::progutils::CmdArgs & inputCommands) {

  std::string directoryName = ".";
  VecStr contains;
  std::string filesOrDirectories = "both";

  bool specific = false;
  bool recursive = false;
  miscSetUp setUp(inputCommands);
  setUp.setUpListAllFiles(directoryName, contains, specific, recursive,
                          filesOrDirectories);
  if (specific) {
    std::cout << "looking for files containing "
              << vectorToString(contains, ", ") << std::endl;
  }
  if (recursive) {
    std::cout << "Looking recursively" << std::endl;
  }
  std::map<std::string, std::pair<std::string, bool>> allFiles = getFiles(
      directoryName, contains, filesOrDirectories, specific, recursive);
  std::map<std::string, std::pair<std::string, bool>>::iterator fileIter;
  // now output the files
  std::cout << std::left;
  uint32_t firstMax = 0;
  uint32_t secondMax = 0;
  for (fileIter = allFiles.begin(); fileIter != allFiles.end(); ++fileIter) {
    if (fileIter->first.length() > firstMax) {
      firstMax = fileIter->first.length();
    }
    if (fileIter->second.first.length() > secondMax) {
      secondMax = fileIter->first.length();
    }
  }
  for (fileIter = allFiles.begin(); fileIter != allFiles.end(); ++fileIter) {
    std::cout << std::setw((int)firstMax) << std::right << fileIter->first
              << " " << std::left << std::setw((int)secondMax)
              << fileIter->second.first << std::endl;
  }

  return 0;
}





int miscRunner::getSharedLines(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	std::string file1 = "";
	std::string file2 = "";
	setUp.setOption(file1, "-file1", "file1");
	setUp.setOption(file2, "-file2", "file2");
	setUp.finishSetUp(std::cout);
	std::unordered_map<std::string, std::string> file1Map;
	std::unordered_map<std::string, std::string> file2Map;
	std::ifstream f1(file1);
	for(std::string line; std::getline(f1, line); ){
		VecStr toks = stringToVector<std::string>(line);
		file1Map[toks[1]] = line;
	}
	std::ifstream f2(file2);
	for(std::string line; std::getline(f2, line); ){
		VecStr toks = stringToVector<std::string>(line);
		file2Map[toks[1]] = line;
	}
	VecStr sharedLines;
	VecStr onlyInOne;
	VecStr onlyInTwo;
	VecStr sharedBuDifferent;
	for(const auto & line : file1Map){
		if(file2Map.find(line.first) == file2Map.end()){
			onlyInOne.emplace_back(line.second);
		}else if(file2Map[line.first] != line.second){
			sharedBuDifferent.emplace_back(line.second + "::::" + file2Map[line.first]);
		}else{
			sharedLines.emplace_back(line.second );
		}
	}
	for(const auto & line : file2Map){
		if(file1Map.find(line.first) == file1Map.end()){
			onlyInTwo.emplace_back(line.second );
		}else{

		}
	}
	std::cout << njh::bashCT::addBGColor(161) << "onlyInOne" << njh::bashCT::reset << std::endl;
	std::cout << njh::bashCT::addColor(161);
	printVector(onlyInOne, "\n");
	std::cout << njh::bashCT::reset;
	std::cout << njh::bashCT::addBGColor(172) << "onlyInTwo" << njh::bashCT::reset << std::endl;
	std::cout << njh::bashCT::addColor(172);
	printVector(onlyInTwo, "\n");
	std::cout << njh::bashCT::reset;
	std::cout << njh::bashCT::addBGColor(41) << "sharedLines" << njh::bashCT::reset << std::endl;
	std::cout << njh::bashCT::addColor(41);
	printVector(sharedLines, "\n");
	std::cout << njh::bashCT::reset;
	std::cout << njh::bashCT::addBGColor(26) << "sharedBuDifferent" << njh::bashCT::reset << std::endl;
	std::cout << njh::bashCT::addColor(26);
	printVector(sharedBuDifferent, "\n");
	std::cout << njh::bashCT::reset;
	return 0;
}

int miscRunner::aminoAcidQuiz(const njh::progutils::CmdArgs & inputCommands) {
  seqSetUp setUp(inputCommands);
  uint32_t limit = 5;
  bool repeats = false;
  std::string previousResults = "";
  setUp.setOption(previousResults, "-previous,-previousResults", "previousResults");
  setUp.setOption(repeats, "-repeats,-repeat", "repeats");
  setUp.setOption(limit, "-limit", "limit");
  setUp.finishSetUp(std::cout);
  uint32_t correct = 0;
  uint32_t incorrect = 0;
  VecStr tests{"Codons", "Single Letter Code", "Triple Letter Code", "Full Name", "Classification"};
  randomGenerator gen;
  std::vector<char> aminos = getAAs(false);
  std::multimap<std::tuple<char,std::string, std::string> , std::pair<uint32_t, bool>> results;
  if(previousResults != ""){
  	table previous(previousResults, "\t", true);
  	for(const auto & row : previous.content_){
  		results.insert({std::make_tuple(row[0][0], row[1], row[2]) , {std::stoi(row[3]),true }});
  	}
  }
  while (correct < limit){
  	std::string given = gen.unifRandSelection(tests);
  	while(given == "Classification"){
  		given = gen.unifRandSelection(tests);
  	}
  	std::string currentTest = given;
  	while(currentTest == given){
  		currentTest = gen.unifRandSelection(tests);
    	while(given == "Codons" && currentTest == "Classification"){
    		currentTest = gen.unifRandSelection(tests);
    	}
  	}
  	char randomAA = gen.unifRandSelection(aminos);
  	std::string give = "";
  	if(given == "Codons"){
  		give = vectorToString(aminoAcidInfo::infos::allInfo.at(randomAA).rnaCodons_, ",");
  	}else if(given == "Single Letter Code"){
  		give = estd::to_string(aminoAcidInfo::infos::allInfo.at(randomAA).letCode_);
  	}else if(given == "Triple Letter Code"){
  		give = aminoAcidInfo::infos::allInfo.at(randomAA).triCode_;
  	}else if(given == "Full Name"){
  		give = aminoAcidInfo::infos::allInfo.at(randomAA).fullName_;
  	}else if(given == "Classification"){
  		give = aminoAcidInfo::infos::allInfo.at(randomAA).classification_;
  	}
  	if(results.find(std::make_tuple(randomAA, given, currentTest)) != results.end() && !repeats){
  		//already asked this question
  		continue;
  	}
  	std::cout << std::endl;
  	std::cout << njh::bashCT::addBGColor(47) << "Correct: " << correct << njh::bashCT::reset << std::endl;
  	std::cout << njh::bashCT::addBGColor(196) << "Incorrect: " << incorrect << njh::bashCT::reset << std::endl;
  	std::cout << "Number of correct guess needed: " << limit - correct << std::endl;
  	std::cout << given << std::endl;
  	std::cout << give << std::endl;
  	std::cout << currentTest << "?" << std::endl;

  	if(currentTest =="Codons"){
  		std::cout << "Give at least one RNA codon, give more by separating by comma" << std::endl;
  		std::string codonAttempt;
  		std::cin >> codonAttempt;
  		VecStr attempt = tokenizeString(codonAttempt, ",");
  		VecStr correctGuesses;
  		VecStr incorrectGuesses;
  		std::string testCodons = vectorToString(aminoAcidInfo::infos::allInfo.at(randomAA).rnaCodons_, ",");
  		for(const auto & tok : attempt){
  			if(testCodons.find(tok) != std::string::npos){
  				correctGuesses.emplace_back(tok);
  			}else{
  				incorrectGuesses.emplace_back(tok);
  			}
  		}
  		if(incorrectGuesses.empty() && !correctGuesses.empty()){
  			std::cout << njh::bashCT::addBGColor(47) << "Correct!" << njh::bashCT::reset << std::endl;
  			std::cout << njh::bashCT::addBGColor(47) << vectorToString(correctGuesses, ",") << njh::bashCT::reset << std::endl;
  			VecStr toksCor = tokenizeString(testCodons, ",");
  			VecStr otherAnswers;
  			for(const auto & cor : toksCor){
  				if(!njh::in(cor, correctGuesses)){
  					otherAnswers.emplace_back(cor);
  				}
  			}
  			if(!otherAnswers.empty()){
  				std::cout << "Other correct answers: " << std::endl;
  				std::cout << njh::bashCT::addBGColor(47) << vectorToString(otherAnswers, ",") << njh::bashCT::reset << std::endl;
  			}
  			++correct;
  			results.insert({std::make_tuple(randomAA, given, currentTest), {1, false}});
  		}else if(incorrectGuesses.empty() && correctGuesses.empty()){
  			std::cout << "no guess recorded" << std::endl;
  		}else{
  			std::cout << njh::bashCT::addBGColor(196) << "incorrect!" << njh::bashCT::reset << std::endl;
  			if(!correctGuesses.empty()){
  				std::cout << "Correct guesses" << std::endl;
  				std::cout << njh::bashCT::addBGColor(47) << vectorToString(correctGuesses, ",") << njh::bashCT::reset << std::endl;
  			}
  			std::cout << njh::bashCT::addBGColor(196) << vectorToString(incorrectGuesses, ",") << njh::bashCT::reset << std::endl;
  			std::cout << "Possible answers were" << std::endl;
  			std::cout << njh::bashCT::addBGColor(47) <<  testCodons << njh::bashCT::reset << std::endl;
  			++incorrect;
  			results.insert({std::make_tuple(randomAA, given, currentTest), {0, false}});
  		}
  	}else if(currentTest =="Single Letter Code"){
  		char guess = ' ';
  		std::cin >> guess;
  		if(guess == ' '){
  			std::cout << "no guess recorded" << std::endl;
  		}else{
  			if(toupper(guess) == aminoAcidInfo::infos::allInfo.at(randomAA).letCode_){
  				std::cout << njh::bashCT::addBGColor(47) << "Correct!" << njh::bashCT::reset << std::endl;
  				++correct;
  				results.insert({std::make_tuple(randomAA, given, currentTest),{1, false}});
  			}else{
  				std::cout << njh::bashCT::addBGColor(196) << "incorrect!" << njh::bashCT::reset << std::endl;
  				std::cout << "answer" << std::endl;
  				std::cout << njh::bashCT::addBGColor(47) << aminoAcidInfo::infos::allInfo.at(randomAA).letCode_ << njh::bashCT::reset << std::endl;
  				++incorrect;
  				results.insert({std::make_tuple(randomAA, given, currentTest), {0, false}});
  			}
  		}
  	}else if(currentTest =="Triple Letter Code"){
  		std::string guess = "";
			std::cin >> guess;
			if(guess == ""){
				std::cout << "no guess recorded" << std::endl;
			}else{
				if(stringToLowerReturn(guess) == stringToLowerReturn(aminoAcidInfo::infos::allInfo.at(randomAA).triCode_)){
					std::cout << njh::bashCT::addBGColor(47) << "Correct!" << njh::bashCT::reset << std::endl;
					++correct;
					results.insert({std::make_tuple(randomAA, given, currentTest),{1, false}});
				}else{
					std::cout << njh::bashCT::addBGColor(196) << "incorrect!" << njh::bashCT::reset << std::endl;
					std::cout << "answer" << std::endl;
					std::cout << njh::bashCT::addBGColor(47) << stringToLowerReturn(aminoAcidInfo::infos::allInfo.at(randomAA).triCode_)
							<< njh::bashCT::reset << std::endl;
					++incorrect;
					results.insert({std::make_tuple(randomAA, given, currentTest), {0, false}});
				}
			}
  	}else if(currentTest =="Full Name"){
  		std::string guess = "";
			std::cin >> guess;
			if(guess == ""){
				std::cout << "no guess recorded" << std::endl;
			}else{
				if(stringToLowerReturn(guess) == stringToLowerReturn(aminoAcidInfo::infos::allInfo.at(randomAA).fullName_)){
					std::cout << njh::bashCT::addBGColor(47) << "Correct!" << njh::bashCT::reset << std::endl;
					++correct;
					results.insert({std::make_tuple(randomAA, given, currentTest),{1, false}});
				}else{
					std::cout << njh::bashCT::addBGColor(196) << "incorrect!" << njh::bashCT::reset << std::endl;
					std::cout << "answer" << std::endl;
					std::cout << njh::bashCT::addBGColor(47) << stringToLowerReturn(aminoAcidInfo::infos::allInfo.at(randomAA).fullName_)
							<< njh::bashCT::reset << std::endl;
					++incorrect;
					results.insert({std::make_tuple(randomAA, given, currentTest), {0, false}});
				}
			}
  	}else if(currentTest =="Classification"){
  		std::cout << "options are : nonpolar, polar, acidic, basic, or stop" << std::endl;
  		std::string guess = "";
			std::cin >> guess;
			if(guess == ""){
				std::cout << "no guess recorded" << std::endl;
			}else{
				if(stringToLowerReturn(guess) == stringToLowerReturn(aminoAcidInfo::infos::allInfo.at(randomAA).classification_)){
					std::cout << njh::bashCT::addBGColor(47) << "Correct!" << njh::bashCT::reset << std::endl;
					++correct;
					results.insert({std::make_tuple(randomAA, given, currentTest),{1, false}});
				}else{
					std::cout << njh::bashCT::addBGColor(196) << "incorrect!" << njh::bashCT::reset << std::endl;
					std::cout << "answer" << std::endl;
					std::cout << njh::bashCT::addBGColor(47) << stringToLowerReturn(aminoAcidInfo::infos::allInfo.at(randomAA).classification_)
							<< njh::bashCT::reset << std::endl;
					++incorrect;
					results.insert({std::make_tuple(randomAA, given, currentTest), {0, false}});
				}
			}
  	}
  }
  std::ofstream aminoAcidQuizResultsFile;
  bool alreadyStarted = bfs::exists("aminoAcidQuizResultsFile.tab.txt");
  aminoAcidQuizResultsFile.open("aminoAcidQuizResultsFile.tab.txt", std::ios::app);
  if(!alreadyStarted){
  	aminoAcidQuizResultsFile << "AminoAcid\tGiven\tTest\tresult" << std::endl;
  }
  for(const auto & result : results ){
  	if(!result.second.second){
  		aminoAcidQuizResultsFile << std::get<0>(result.first)
  		  			<< "\t" << std::get<1>(result.first)
  		  			<< "\t" << std::get<2>(result.first)
  		  			<< "\t" << result.second.first << std::endl;
  	}
  }
  return 0;
}









}  // namespace njh
