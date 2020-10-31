//
// Created by Andrew Bailey on 6/7/20.
//

#ifndef EMBED_FAST5_SRC_AMBIGMODEL_HPP_
#define EMBED_FAST5_SRC_AMBIGMODEL_HPP_

#include "EmbedUtils.hpp"

using namespace embed_utils;
using namespace std;

class AmbigModel {
 public:
  std::map<string, string> ambig_map;
  std::map<char, std::set<char>> canonical_counterpart;
  std::map<char, std::set<char>> mod_counterparts;

  std::set<char> canonical{'A', 'C', 'G', 'T'};

  AmbigModel(string ambig_model) {
    load(ambig_model);
  }
  ~AmbigModel() = default;
  AmbigModel() = default;

  void load(string ambig_model){
    std::map<string, string> ambig_hash;
    char ambig_bases[10];
    char nucleotides[10];
    char line[100];
    FILE *infile = fopen(ambig_model.c_str(), "r");
    throw_assert(infile, "Couldn't open " + string(ambig_model) + " for reading\n")
    int i = 0;
    while(i < 300 && fgets(line, sizeof(line), infile) != nullptr){
      sscanf(line, "%s\t%s", ambig_bases, nucleotides);
      ambig_map.insert(std::pair<string, string>(ambig_bases, nucleotides));
      this->parse_nucleotides_from_model(nucleotides);
      i++;
    }
  }

  void parse_nucleotides_from_model(const string& nucleotides){
    string canonical_bases;
    string mod;
    for (char i: nucleotides) {
      if (is_canonical(i)) {
        canonical_bases += i;
      } else {
        mod += i;
      }
    }
    if (canonical_bases.size() > 0){
      for (char &i: mod){
        canonical_counterpart[i].insert(canonical_bases.begin(), canonical_bases.end());
      }
    }
    for (char &i: canonical_bases){
      mod_counterparts[i].insert(mod.begin(), mod.end());
    }
  }

  std::set<string> get_canonical_kmers(const string& kmer){
    vector<string> canonical_kmers{""};
    for (auto &i: kmer){
      if (is_canonical(i)){
        for (auto &k: canonical_kmers){
          k += i;
        }
      } else {
        throw_assert(in_canonical_counterpart(i), "Character("+ string(1, i) +")in kmer not in ambig model");
        std::set<char> canonical_chars = canonical_counterpart[i];
        vector<string> intermediate_kmers;
        for (auto &c: canonical_chars){
          vector<string> tmp_kmers = canonical_kmers;
          for (auto &k: tmp_kmers){
            k += c;
          }
          intermediate_kmers.insert(intermediate_kmers.end(), tmp_kmers.begin(),tmp_kmers.end());
        }
        canonical_kmers = intermediate_kmers;
      }
    }
    std::set<string> s;
    for (auto x : canonical_kmers) {
      s.insert(x);
    }
    return s;
  }

  std::set<string> get_canonical_kmers(const std::set<string>& kmers){
    std::set<string> s;
    for (auto &kmer: kmers){
      std::set<string> a = get_canonical_kmers(kmer);
      s.insert(a.begin(), a.end());
    }
    return s;
  }

  bool is_canonical(const char& nuc){
    return canonical.find(nuc) != canonical.end();
  }

  bool is_canonical(const string& kmer){
    for (auto &nuc: kmer){
      if (!is_canonical(nuc)){
        return false;
      }
    }
    return true;
  }

  bool is_canonical(const std::set<string>& canonical){
    for (auto &i: canonical){
      if (!is_canonical(i)){
        return false;
      }
    }
    return true;
  }

  bool in_canonical_counterpart(const char& nuc){
    return canonical_counterpart.find(nuc) != canonical_counterpart.end();
  }
  bool in_mod_counterpart(const char& nuc){
    return mod_counterparts.find(nuc) != mod_counterparts.end();
  }

};

#endif //EMBED_FAST5_SRC_AMBIGMODEL_HPP_
