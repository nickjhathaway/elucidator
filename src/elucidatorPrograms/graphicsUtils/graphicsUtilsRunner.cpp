//
//  graphicsUtils.cpp
//  sequenceTools
//
//  Created by Nicholas Hathaway on 05/30/2013.
//  Copyright (c) 2013 Nicholas Hathaway. All rights reserved.
//

#include "graphicsUtilsRunner.hpp"


namespace njhseq {


graphicsUtilsRunner::graphicsUtilsRunner()
    : njh::progutils::ProgramRunner({addFunc("printAnsiColors", printAnsiColors, false),
										 addFunc("colorInfo", colorInfo, false),
                     addFunc("getColors", getColors, false),
                     addFunc("multipleColorsInfo", multipleColorsInfo, false)},
                    "graphicsUtils") {}


int graphicsUtilsRunner::printAnsiColors(const njh::progutils::CmdArgs & inputCommands) {
	seqSetUp setUp(inputCommands);
	setUp.finishSetUp(std::cout);
	std::multimap<std::tuple<double, double, double>,
			std::pair<njh::color, uint32_t>> byHueColor;
	auto ansiColorCodeToColor = njh::colorspace::createAnsiColorCodeToColor();
	for (const auto & ansiCol : ansiColorCodeToColor) {
		byHueColor.insert(
				std::make_pair(
						std::tuple<double, double, double>(ansiCol.second.lum_,
								ansiCol.second.lSat_, ansiCol.second.hue_),
						std::make_pair(ansiCol.second, ansiCol.first)));
	}
  njh::color lastColor = byHueColor.rbegin()->second.first;
  table colorInfoOut(getColorInfo(lastColor, njh::color(), true));
  for(const auto & h : byHueColor){
  	colorInfoOut.content_.emplace_back(getColorInfo(h.second.first, lastColor, false, true));
  	lastColor = h.second.first;
  }
  colorInfoOut.outPutContentOrganized(std::cout);
  return 0;
}

int graphicsUtilsRunner::colorInfo(const njh::progutils::CmdArgs & inputCommands) {
  seqSetUp setUp(inputCommands);
  setUp.finishSetUp(std::cout);
  std::string colorMode = "";
  while (colorMode!= "rgb" && colorMode != "hsl" && colorMode != "hex"){
  	std::cout << "color by rgb, hsl, or hex" << std::endl;
  	std::cin >> colorMode;
  	stringToLower(colorMode);
  }
  njh::color outColor;
  if(colorMode == "rgb"){
  	uint32_t maxColor = 0;
  	while(maxColor!= 1 && maxColor!= 255){
  		std::cout << "Enter maxColor, options are 1 or 255" << std::endl;
  		std::cin >> maxColor;
  	}
  	if(maxColor == 1){
    	double redValue = 1.01;
    	double blueValue = 1.01;
    	double greenValue = 1.01;
    	while(redValue > 1 || redValue < 0){
    		std::cout << "Enter red value, should be between 0 and 1" << std::endl;
    		std::cin >> redValue;
    	}
    	while(greenValue > 1 || greenValue < 0){
    		std::cout << "Enter green value, should be between 0 and 1" << std::endl;
    		std::cin >> greenValue;
    	}
    	while(blueValue > 1 || blueValue < 0){
    		std::cout << "Enter blue value, should be between 0 and 1" << std::endl;
    		std::cin >> blueValue;
    	}
    	outColor = njh::color(redValue, greenValue, blueValue, maxColor);
  	}else{
    	uint32_t redValue = 256;
    	uint32_t blueValue = 256;
    	uint32_t greenValue = 256;
    	while(redValue > 255){
    		std::cout << "Enter red value, should be between 0 and 255" << std::endl;
    		std::cin >> redValue;
    	}
    	while(greenValue > 255){
    		std::cout << "Enter green value, should be between 0 and 255" << std::endl;
    		std::cin >> greenValue;
    	}
    	while(blueValue > 255){
    		std::cout << "Enter blue value, should be between 0 and 255" << std::endl;
    		std::cin >> blueValue;
    	}
    	outColor = njh::color(redValue, greenValue, blueValue, maxColor);
  	}
  }else if (colorMode == "hsl"){
  	double hueValue = 361;
  	double satValue = 1.01;
  	double lumValue = 1.01;
  	while(hueValue > 360 || hueValue < 0){
  		std::cout << "Enter hue value, should be between 0 and 360" << std::endl;
  		std::cin >> hueValue;
  	}
  	while(satValue > 1 || satValue < 0){
  		std::cout << "Enter saturation value, should be between 0 and 1" << std::endl;
  		std::cin >> satValue;
  	}
  	while(lumValue > 1 || lumValue < 0){
  		std::cout << "Enter luminance value, should be between 0 and 1" << std::endl;
  		std::cin >> lumValue;
  	}
  	outColor.setColorByHueSatLum(hueValue, satValue, lumValue);
  } else if (colorMode == "hex"){
  	std::string hexString = "";
  	while(hexString.length() <6 || hexString.length() > 9){
  		std::cout << "Enter hex string for color, can start with # or not" << std::endl;
  		std::cout << "Needs to be either 6 letters or 8 letters (if including alpha)\n";
  		std::cout << "Eg. for black #FFFFFF (or FFFFFF) or black with alpha #FFFFFFFF (or FFFFFFFF)\n";
  		std::cin >> hexString;
  	}
  	//std::cout << hexString << std::endl;
  	outColor.setColorByHex(hexString);
  }
  table colorInfoOut(getColorInfo(outColor, true));
  colorInfoOut.content_.emplace_back(getColorInfo(outColor));
  colorInfoOut.outPutContentOrganized(std::cout);
  return 0;
}

int graphicsUtilsRunner::getColors(const njh::progutils::CmdArgs & inputCommands) {
  seqSetUp setUp(inputCommands);
  bool grayScale = false;
  double satStart = 1.00;
  double satStop = 1.00;
  double lumStart = 0.6;
  double lumStop = 0.6;
  double alpha = 1.00;
  bool alphaSet = false;
  double hueStart = 0.0;
  double hueStop = 360.0;

  uint32_t colorNumber = 10;
  bool rOutPut = false;
  std::string pdfFilename = "";
  setUp.setOption(pdfFilename, "-pdfFilename,-pdf,-out", "pdfFilename");
  setUp.setOption(rOutPut, "-r,-rOutput", "rOutPut");
  setUp.setOption(colorNumber, "-number,-num,-n", "colorNumber");
  setUp.setOption(grayScale, "-gray,-grayScale", "grayScale");
  setUp.setOption(satStart, "-sat,-satStart", "satStart");
  if(!setUp.setOption(satStop, "-satStop", "satStop")){
  	satStop = satStart;
  }
  setUp.setOption(lumStart, "-lum,-lumStart", "lumStart");
  if(!setUp.setOption(lumStop, "-lumStop", "lumStop")){
  	lumStop = lumStart;
  }
  setUp.setOption(hueStart, "-hueStart", "hueStart");
  setUp.setOption(hueStop, "-hueStop", "hueStop");
  if(setUp.setOption(alpha, "-alpha", "alpha")){
  	alphaSet = true;
  }

  setUp.finishSetUp(std::cout);

  if(hueStart == 0.0 && hueStop == 360.0){
  	double hueStep = hueStop / colorNumber;
  	hueStop -= hueStep;
  }
  std::vector<njh::color> outColors;
  outColors = njh::getColsBetweenInc(hueStart, hueStop, lumStart, lumStop, satStart, satStop, colorNumber);
  if(alphaSet){
  	double maxAlpha = 1;
  	if(alpha > 1){
  		maxAlpha = 255;
  	}
  	for_each(outColors, [&alpha, maxAlpha](njh::color & col){ col.setAlpha(alpha, maxAlpha);});
  }
  if(grayScale){
  	for(auto pos : iter::range(outColors.size())){
  		outColors[pos] = outColors[pos].getGreyScaleColor();
  	}
  }

  if(rOutPut){
  	printColorsForR(outColors, std::cout);
  }else{
    //table colorInfoOut(VecStr{"hexStr", "hue", "sat", "lum", "r,g,b", "lastDiff" ,"closestWord", "closestAnsi"});
  	table colorInfoOut(getColorInfo(njh::color(), njh::color(), true, true));
    njh::color lastColor = outColors.back();
    for(const auto & col : outColors){
    	auto closeAnsi = getClosetAnsiColor(col);
    	colorInfoOut.content_.emplace_back(getColorInfo(col, lastColor, false, true));
    	lastColor = col;
    }
    colorInfoOut.outPutContentOrganized(std::cout);
  }

  return 0;
}

int graphicsUtilsRunner::multipleColorsInfo(const njh::progutils::CmdArgs & inputCommands) {
	std::string inColorsStr = "#F7FBFF,#DEEBF7,#C6DBEF,#9ECAE1,#6BAED6,#4292C6,#2171B5,#08519C,#08306B";
  graphicsUtilsSetUp setUp(inputCommands);
  setUp.setOption(inColorsStr, "-in,-inColorsStr,-colors", "inColorsStr");
  setUp.finishSetUp(std::cout);
  VecStr inColors = tokenizeString(inColorsStr, ",");
  table colorInfo(getColorInfo(njh::color(), njh::color(), true, true));
  auto outColors = njh::multipleHexToColor(inColors);
  njh::color lastColor = outColors.back();
  for(const auto & col : outColors){
  	colorInfo.content_.emplace_back(getColorInfo(col, lastColor, false, true));
  	lastColor = col;
  }
  colorInfo.outPutContentOrganized(std::cout);
  return 0;
}




}  // namespace njh
