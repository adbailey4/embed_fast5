//
// Created by Andrew Bailey on 5/2/20.
// Most of this code came from Ryan Lorig-Roach. https://github.com/rlorigro/runlength_analysis_cpp
//

//#include "BinaryEventWriter.hpp"

//BinaryRunnieWriter::BinaryRunnieWriter(path file_path) {
//  this->sequence_file_path = file_path;
//
//  // Ensure that the output directory exists
//  create_directories(this->sequence_file_path.parent_path());
//
//  this->sequence_file = ofstream(this->sequence_file_path,ofstream::binary);
//
//  if (not this->sequence_file.is_open()){
//    throw runtime_error("ERROR: could not open file " + file_path.string());
//  }
//}
//
//
//void BinaryRunnieWriter::write_sequence_block(RunnieSequenceElement& sequence){
//  // Write the sequence to the file
//  write_string_to_binary(this->sequence_file, sequence.sequence);
//}
//
//
//void BinaryRunnieWriter::write_scales_block(RunnieSequenceElement& sequence){
//  // Write the encodings to the file
//  for (auto& length: sequence.scales){
//    write_value_to_binary(this->sequence_file, length);
//  }
//}
//
//
//void BinaryRunnieWriter::write_shapes_block(RunnieSequenceElement& sequence){
//  // Write the encodings to the file
//  for (auto& length: sequence.shapes){
//    write_value_to_binary(this->sequence_file, length);
//  }
//}
//
//
//void BinaryRunnieWriter::write_sequence(RunnieSequenceElement& sequence){
//  if (sequence.sequence.empty()){
//    throw runtime_error("ERROR: empty sequence provided to BinaryRunnieWriter: " + sequence.name);
//  }
//
//  RunlengthIndex index;
//
//  // Add sequence start position to index
//  index.sequence_byte_index = this->sequence_file.tellp();
//
//  // Store the length of this sequence
//  index.sequence_length = sequence.sequence.size();
//
//  // Store the name of this sequence
//  index.name = sequence.name;
//
//  this->write_sequence_block(sequence);
//  this->write_scales_block(sequence);
//  this->write_shapes_block(sequence);
//
//  // Append index object to vector
//  this->indexes.push_back(index);
//}
//
//
//void BinaryRunnieWriter::write_index(RunlengthIndex& index){
//  // Where is the sequence
//  write_value_to_binary(this->sequence_file, index.sequence_byte_index);
//
//  // How long is the sequence
//  write_value_to_binary(this->sequence_file, index.sequence_length);
//
//  // How long is the name of the sequence
//  write_value_to_binary(this->sequence_file, index.name.size());
//
//  // What is the name
//  write_string_to_binary(this->sequence_file, index.name);
//}
//
//
//void BinaryRunnieWriter::write_indexes(){
//  // Store the current file byte index so the beginning of the INDEX table can be located later
//  uint64_t indexes_start_position = this->sequence_file.tellp();
//
//  // Iterate all the indexes, write them to the file
//  for (auto& index: this->indexes){
//    write_index(index);
//  }
//
//  // Store the current file byte index so the beginning of the CHANNEL table can be located later
//  uint64_t channel_metadata_start_position = this->sequence_file.tellp();
//
//  // Write channel metadata
//  write_value_to_binary(this->sequence_file, BinaryRunnieWriter::n_channels);
//  write_value_to_binary(this->sequence_file, BinaryRunnieWriter::channel_sizes[BinaryRunnieWriter::LENGTH]);
//
//  // Write the pointer to the beginning of the index table
//  write_value_to_binary(this->sequence_file, indexes_start_position);
//
//  // Write the pointer to the beginning of the channels table
//  write_value_to_binary(this->sequence_file, channel_metadata_start_position);
//}
