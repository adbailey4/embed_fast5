//
// Created by Andrew Bailey on 2019-07-26.
//

#ifndef EMBED_FAST5_SRC_FOLDERHANDLER_HPP_
#define EMBED_FAST5_SRC_FOLDERHANDLER_HPP_

// Embed
#include "EmbedUtils.hpp"

//Boost Libraries
#include <boost/filesystem.hpp>
#include <boost/coroutine2/all.hpp>

// Standard Libraries
#include <string>

using namespace boost::filesystem;
using namespace boost::coroutines2;
using namespace std;
using namespace embed_utils;


//class FolderHandler{
// public:
//  FolderHandler(path input_dir, uint64_t batch_size, string ext);
//  ~FolderHandler();
//  path input_dir;
//  uint64_t batch_size;
//  string ext;
//  vector<string> get_batch();
////  boost::coroutines2::detail::pull_coroutine<path> generator;
// private:
//  dir_coro::pull_type* generator{};
//};

#endif //EMBED_FAST5_SRC_FOLDERHANDLER_HPP_
