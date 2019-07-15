//
// Created by Andrew Bailey on 2019-07-15.
//

#include "AlignmentFile.hpp"
#include "EmbedUtils.hpp"

using namespace embed_utils;

AlignmentFile::~AlignmentFile()
= default;

AlignmentFile::AlignmentFile(const string& input_reads_filename){
  this->file_path = input_reads_filename;
  this->get_strand();
  this->get_k();
}

int64_t AlignmentFile::get_k(){
  std::ifstream in_file(this->file_path.c_str());
  if (in_file.good()) {
    // read the file
    std::string line;
    getline(in_file, line);
    std::vector<std::string> fields = split(line, '\t');
    this->k = fields[2].length();
  }
  in_file.close();
  return this->k;
}

void AlignmentFile::filter(PositionsFile* pf, boost::filesystem::path& output_file, string bases) {
  std::ofstream out_file;
  out_file.open(output_file.string());
  std::ifstream in_file(this->file_path.c_str());

  if (in_file.good()) {
    // read the file
    std::string line;
    while (getline(in_file, line)) {
      std::vector<std::string> fields = split(line, '\t');
      string contig = fields[0];
      int64_t reference_index = convert_to_int(fields[1]);
//            string reference_kmer = fields[2];
//            string read_file = fields[3];
      string read_strand = fields[4];
//            string event_index = fields[5];
//            string event_mean = fields[6];
//            string event_noise = fields[7];
//            string event_duration = fields[8];
//            string aligned_kmer = fields[9];
//            string scaled_mean_current = fields[10];
//            string scaled_noise = fields[11];
      string posterior_probability = fields[12];
      string descaled_event_mean = fields[13];
//            string ont_model_mean = fields[14];
      string path_kmer = fields[15];
      string contig_strand = contig+this->strand;
      if (pf->is_in(contig_strand, reference_index) || !are_characters_in_string(bases, path_kmer)) {
        out_file << path_kmer << '\t' << read_strand << '\t' << descaled_event_mean << '\t' <<  posterior_probability << '\n';
      }
    }
  }
  out_file.close();
  in_file.close();
}



string AlignmentFile::get_strand(){
  std::vector<std::string> fields = split(this->file_path, '.');
  if (fields.end()[-2] == "backward"){
    this->strand = "-";
  } else if (fields.end()[-2] == "forward") {
    this->strand = "+";
  } else {
    fprintf(stderr, "error: could not infer strand from  %s\n", this->file_path.c_str());
    fprintf(stderr, "Please check input file is full alignment file from signalalign\n");
    exit(EXIT_FAILURE);

  }
  return this->strand;
}