//
// Created by Nicholas Hathaway on 8/20/25.
//

#include "pmo_utils_runner.hpp"


#include <njhseq/IO/SeqIO/SeqIO.hpp>
#include <njhseq/objects/BioDataObject/pmo.h>
#include <nlohmann/json.hpp>

namespace njhseq {
PMOUtilsRunner::PMOUtilsRunner()
        : njh::progutils::ProgramRunner(
        {
        	addFunc("read_pmo", read_pmo, false),
        	addFunc("get_overlap_between_panels_in_pmos", get_overlap_between_panels_in_pmos, false),
        	addFunc("add_protein_variant_info_to_pmo", add_protein_variant_info_to_pmo, false),

        },//
        "PMOUtils") {}


int PMOUtilsRunner::read_pmo(const njh::progutils::CmdArgs & inputCommands) {
	bfs::path pmo_fnp;
	OutOptions out_options;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(pmo_fnp, "--pmo", "pmo file to read", true);
	setUp.processWritingOptions(out_options);
	setUp.finishSetUp(std::cout);

	InputStream in(pmo_fnp);
	auto pmo_obj= pmo::PortableMicrohaplotypeObject::from_json(nlohmann::json::parse(in));
	pmo_obj.validate();
	OutputStream out(out_options);
	out << pmo_obj.to_json();

	return 0;
}


int PMOUtilsRunner::get_overlap_between_panels_in_pmos(const njh::progutils::CmdArgs & inputCommands) {
	std::vector<bfs::path> pmo_fnps;
	bool strands_must_match = false;
	OutOptions out_options;
	seqSetUp setUp(inputCommands);
	setUp.processVerbose();
	setUp.processDebug();
	setUp.setOption(pmo_fnps, "--pmos", "A list of file name paths to read in", true);
	setUp.setOption(strands_must_match, "--strands_must_match", "strands must match");

	setUp.processWritingOptions(out_options);
	setUp.finishSetUp(std::cout);
	std::vector<pmo::PortableMicrohaplotypeObject> pmo_objs;
	OutputStream out(out_options);

	for (const auto & pmo_fnp : pmo_fnps) {
		InputStream in(pmo_fnp);
		pmo_objs.emplace_back(pmo::PortableMicrohaplotypeObject::from_json(nlohmann::json::parse(in)));
		pmo_objs.back().validate();
	}
	out << "panel1";
	out << "\t" << "panel1_chrom\tpanel1_start\tpanel1_end\tpanel1_name\tpanel1_len\tpanel1_strand";

	out << "\t" << "panel2";
	out << "\t" << "panel2_chrom\tpanel2_start\tpanel2_end\tpanel2_name\tpanel2_len\tpanel2_strand";

	out << "\t" << "overlap_chrom\toverlap_start\toverlap_end\toverlap_name\toverlap_len";

	out << "\t" << "panel1_coverage" << "\t" << "panel2_coverage" << std::endl;
	std::unordered_map<std::string, std::vector<std::shared_ptr<Bed6RecordCore>>> bedsByPanel;
	for (const auto & pmo_obj : pmo_objs) {
		for (const auto & pmo_panel : pmo_obj.panel_info_) {
			if (njh::notIn(pmo_panel.panel_name_, bedsByPanel)) {
				std::unordered_map<uint32_t, VecStr> targets_to_reactions;
				for (const auto & reaction : pmo_panel.reactions_) {
					for (const auto & target : reaction.panel_targets_) {
						targets_to_reactions[target].emplace_back(reaction.reaction_name_);
					}
				}
				std::vector<std::shared_ptr<Bed6RecordCore>> current_beds;
				for (const auto & target : targets_to_reactions) {
					if (!pmo_obj.target_info_[target.first].insert_location_.has_value()) {
						std::stringstream ss;
						ss << __PRETTY_FUNCTION__ << " " << __FILE__ << " " << __LINE__ << ", error " << "target: " << pmo_obj.target_info_[target.first].target_name_ << " doesn't have isnert location loaded" << "\n";
						throw std::runtime_error{ss.str()};
					}
					current_beds.emplace_back(std::make_shared<Bed6RecordCore>(pmo_obj.target_info_[target.first].insert_location_.value().chrom_,
						pmo_obj.target_info_[target.first].insert_location_.value().start_,
						pmo_obj.target_info_[target.first].insert_location_.value().end_,
						pmo_obj.target_info_[target.first].target_name_,
						pmo_obj.target_info_[target.first].insert_location_.value().end_ - pmo_obj.target_info_[target.first].insert_location_.value().start_,
						pmo_obj.target_info_[target.first].insert_location_.value().strand_.has_value() ? pmo_obj.target_info_[target.first].insert_location_.value().strand_.value().front() : '+'
						));
					auto reactions_str = njh::pasteAsStr("[reactions=",  njh::conToStr(target.second, ","), "]");
					current_beds.back()->extraFields_.emplace_back(reactions_str);
				}
				bedsByPanel[pmo_panel.panel_name_] = current_beds;
			}
		}
	}

	std::vector<std::string> panelNames = getVectorOfMapKeys(bedsByPanel);

	njh::sort(panelNames);
	for(const auto & currentPanelPos : iter::range(panelNames.size())){
		const auto & current_panel_name = panelNames[currentPanelPos];
		const auto & current_panel_regions = bedsByPanel.at(current_panel_name);
		for(const auto & otherPanelPos : iter::range(currentPanelPos + 1, panelNames.size())){
			const auto & other_panel_regions = bedsByPanel[panelNames[otherPanelPos]];
			for(const auto current_bed_pos : iter::range(current_panel_regions.size())){
				const auto & current_bed = current_panel_regions[current_bed_pos];
				std::vector<std::shared_ptr<Bed6RecordCore>> overlapping_regions_other_panel;
				for(const auto & other_panel_bed : other_panel_regions){
					if(current_bed->overlaps(*other_panel_bed, 1) && (!strands_must_match || current_bed->strand_ == other_panel_bed->strand_)){
						overlapping_regions_other_panel.emplace_back(other_panel_bed);
					}
				}
				for(const auto & overlapping_region_other_panel : overlapping_regions_other_panel){
					out << current_panel_name;
					out << "\t" << current_bed->toDelimStr();
					out << "\t" << panelNames[otherPanelPos];
					out << "\t" << overlapping_region_other_panel->toDelimStr();
					uint32_t overlapStart = std::max(overlapping_region_other_panel->chromStart_, current_bed->chromStart_);
					uint32_t overlapStop = std::min(overlapping_region_other_panel->chromEnd_, current_bed->chromEnd_);
					uint32_t overlapLen = overlapStop - overlapStart;
					out << "\t" << overlapping_region_other_panel->chrom_ << "\t" << overlapStart << "\t" << overlapStop << "\t"
							<< njh::pasteAsStr(overlapping_region_other_panel->chrom_, "-", overlapStart, "-", overlapStop) << "\t" << overlapLen;
					out << "\t" << static_cast<double>(overlapLen)/current_bed->length() << "\t" << static_cast<double>(overlapLen)/overlapping_region_other_panel->length() << std::endl;
				}
			}
		}
	}

	return 0;
}


} //namespace njhseq


