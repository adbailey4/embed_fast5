//
// Created by Andrew Bailey on 2019-07-26.
//


// Embed
#include "FolderHandler.hpp"
#include "EmbedUtils.hpp"

// Boost libraries.
#include <boost/filesystem.hpp>

// Standard library.
#include <string>

using namespace boost::filesystem;
using namespace boost::coroutines2;
using namespace std;
using namespace embed_utils;


//FolderHandler::FolderHandler(path input_dir, uint64_t batch_size, string ext) :
//    input_dir(std::move(input_dir)), batch_size(batch_size), ext(std::move(ext))
//{
//  dir_coro::pull_type wp = list_files_in_dir(this->input_dir, this->ext);
//
////  shared_ptr<dir_coro::pull_type> sp = wp.lock();
////
////
////  generator = &tmp;
//}
//
//FolderHandler::~FolderHandler(){
//  delete(generator);
//}
//
//vector<string> FolderHandler::get_batch(){
//  vector<string> batch;
//  batch.resize(this->batch_size);
//  path something = generator->get();
//  cout << something << "\n";
//
//}
