//
// Created by Andrew Bailey on 5/3/20.
// Most of this code came from Ryan Lorig-Roach. https://github.com/rlorigro/runlength_analysis_cpp
//

//#include "BinaryEventReader.hpp"



//size_t BinaryEventReader::get_read_count(){
//  return this->indexes.size();
//}
//
//
//const string& BinaryEventReader::get_read_name(uint64_t read_number){
//  return this->indexes.at(read_number).name;
//}
//
//
//uint64_t BinaryEventReader::get_length(uint64_t read_number){
//  return this->indexes.at(read_number).sequence_length;
//}
//
//RunnieSequenceElement BinaryEventReader::generate_sequence_container(){
//  return RunnieSequenceElement();
//}
//
//
//const string& BinaryEventReader::get_file_name(){
//  return this->sequence_file_path;
//}
//
//
//void BinaryEventReader::get_sequence(RunnieSequenceElement& sequence, uint64_t read_number){
//  sequence = {};
//  off_t byte_index = off_t(this->indexes.at(read_number).sequence_byte_index);
//  pread_string_from_binary(this->sequence_file_descriptor, sequence.sequence, this->indexes.at(read_number).sequence_length, byte_index);
//  pread_vector_from_binary(this->sequence_file_descriptor, sequence.scales, this->indexes.at(read_number).sequence_length, byte_index);
//  pread_vector_from_binary(this->sequence_file_descriptor, sequence.shapes, this->indexes.at(read_number).sequence_length, byte_index);
//}
//
//
//void BinaryEventReader::get_sequence(RunnieSequenceElement& sequence, string& read_name){
//  sequence = {};
//  uint64_t read_number;
//  off_t byte_index;
//
//  try {
//    read_number = this->index_map.at(read_name);
//  }catch(std::out_of_range) {
//    cerr << "\nERROR: " << read_name << " not found in index for file " << this->sequence_file_path << '\n';
//    exit(1);
//  }
//
//  byte_index = off_t(this->indexes.at(read_number).sequence_byte_index);
//
//  pread_string_from_binary(this->sequence_file_descriptor, sequence.sequence, this->indexes.at(read_number).sequence_length, byte_index);
//  pread_vector_from_binary(this->sequence_file_descriptor, sequence.scales, this->indexes.at(read_number).sequence_length, byte_index);
//  pread_vector_from_binary(this->sequence_file_descriptor, sequence.shapes, this->indexes.at(read_number).sequence_length, byte_index);
//}



