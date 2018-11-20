#pragma once
/*
 * GzSimpleBinFile.hpp
 *
 *  Created on: Oct 18, 2016
 *      Author: nick
 */

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

#include "elucidator/common.h"


namespace njhseq {

/**@brief Template wrapper for reading binary data from a gz file where data is stored in chunks of the same size (SIZE) of the same type of data(T)
 *
 */
template<typename T, uint32_t SIZE>
class GzSimpleBinFile {
	gzFile gzInFile_; /**< gz file handler for reading */
	const uint32_t numBytes_ = SIZE * sizeof(T); /**< the number of bytes that need to be read in for each chunck */
public:
	/**@brief a small class to wrap an array of the type and size that need to be read, just a convience
	 *
	 */
	class GzSimpleBinFileDCon {
	public:
		std::array<T, SIZE> data_; /**< data */
		/**@brief fill the data_ with a certain value
		 *
		 * @param value the "empty" value defaults to max of the type T
		 */
		void fillWithEmpty(T value = std::numeric_limits<T>::max()) {
			data_.fill(value);
		}
	};

	/**@brief construct with file name to be read from
	 *
	 * @param fnp
	 */
	GzSimpleBinFile(const bfs::path& fnp): fnp_(fnp){
		if(!bfs::exists(fnp)){
			std::stringstream ss;
			ss<< __PRETTY_FUNCTION__ << ": error, file " << fnp << " doesn't exist" << "\n";
			throw std::runtime_error{ss.str()};
		}
		gzInFile_ = gzopen(fnp_.c_str(), "r");
#if ZLIB_VERNUM >= 0x1280
		gzbuffer(gzInFile_, 128 * 1024);
#endif
	}

	const bfs::path fnp_; /**< file being read from*/

	/**@brief seek in the gzInFile_ handler
	 *
	 * @param pos the file position to seek to
	 * @return whether the seek was successful
	 */
	bool seek(uint64_t pos){
		auto success = gzseek(gzInFile_, pos, SEEK_SET);
		if(-1 == success){
			return false;
		}
		return true;
	}
	/**@brief read in and store data chunk in a GzSimpleBinFileDCon class
	 *
	 * @param d
	 * @return
	 */
	bool read(GzSimpleBinFileDCon & d) {
		//d.clear();
		int bytes_read = gzread(gzInFile_, d.data_.data(), numBytes_);
		if (0 == bytes_read && gzeof(gzInFile_)) {
			d.fillWithEmpty();
			return false;
		}
		if (bytes_read != static_cast<int64_t>(numBytes_)) {
			std::stringstream ss;
			ss << __PRETTY_FUNCTION__ << ": Error in reading: " << fnp_ << std::endl;
			ss << "Read in " << bytes_read << " when expected " << numBytes_
					<< std::endl;
			throw std::runtime_error { ss.str() };
		}
		return true;
	}

	/**@brief factor for data container that can be used to read repeatedly from the gz file
	 *
	 * @return a GzSimpleBinFileDCon with the size and type of the currently template class
	 */
	GzSimpleBinFileDCon genDataContainer() {
		return GzSimpleBinFileDCon { };
	}

	/**@brief close the gzInFile_ on exit
	 *
	 */
	~GzSimpleBinFile(){
		gzclose(gzInFile_);
	}
};

}  // namespace njhseq


