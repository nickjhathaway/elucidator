//
// Created by Nicholas Hathaway on 8/7/25.
//



#include "seqUtilsModRunner.hpp"
#include <njhseq/objects/counters/DNABaseCounter.hpp>
#include <njhseq/IO/SeqIO/SeqIO.hpp>
#include <njhseq/objects/seqObjects/seqKmers.h>
#include <njhseq/readVectorManipulation/readVectorHelpers/readVecSorter.hpp>


namespace njhseq {

int seqUtilsModRunner::generateAllPossibleNucleotidePossibleFromProtein(const njh::progutils::CmdArgs & inputCommands) {
	VecStr input_proteins;
	OutOptions outOpts;
	seqSetUp setUp(inputCommands);
	setUp.description_ = "Generate all possible nucleotide sequences possible for input protein(s) using standard codon libary";
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(input_proteins, "--input_proteins", "the input proteins to convert");
	setUp.processWritingOptions(outOpts);
	setUp.finishSetUp(std::cout);

	OutputStream out(outOpts);

	// convert to upper and check for nonstandard amino acids
	std::set<char> non_standard_aas;
	for (auto &protein: input_proteins) {
		njh::strToUpper(protein);
		for (const auto aa : protein) {
			if (njh::notIn(aa, aminoAcidInfo::infos::allInfo)) {
				non_standard_aas.emplace(aa);
			}
		}
	}
	if (!non_standard_aas.empty()) {
		std::stringstream ss;
		ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << " found the following non-standard amino acids " << "\n";
		ss << njh::conToStr(non_standard_aas, ",") << "\n";
		throw std::runtime_error{ss.str()};
	}

	for (const auto &protein: input_proteins) {
		VecStr output;
		for (const auto aa : protein) {
			if (output.empty()) {
				for (const auto & dna_codon : aminoAcidInfo::infos::allInfo.at(aa).dnaCodons_) {
					output.emplace_back(dna_codon);
				}
			} else {
				VecStr new_output;
				for (auto current_out : output) {
					for (const auto & dna_codon : aminoAcidInfo::infos::allInfo.at(aa).dnaCodons_) {
						new_output.emplace_back(current_out + dna_codon);
					}
				}
				output.swap(new_output);
			}
		}
		out << njh::conToStr(output, "\n") << std::endl;
	}

	return 0;
}


} //namespace njhseq



