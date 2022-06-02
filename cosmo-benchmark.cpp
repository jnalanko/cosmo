#include "SeqIO/SeqIO.hh"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <libgen.h> // basename
#include <sys/mman.h> // mlockall

#include "tclap/CmdLine.h"

#include <sdsl/bit_vectors.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/wt_algorithm.hpp>

#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/iterator/function_input_iterator.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem.hpp>
#include <boost/optional.hpp>
//#include <boost_iterator/zip_iterator.hpp>

#include "io.hpp"
#include "debruijn_graph.hpp"
#include "debruijn_graph_shifted.hpp"
#include "debruijn_hypergraph.hpp"
#include "algorithm.hpp"
#include "wt_algorithm.hpp"


#include <string>
#include <sstream>

using namespace std;
using namespace sdsl;
using namespace std::chrono;

long long cur_time_millis(){
  return (std::chrono::duration_cast< milliseconds >(high_resolution_clock::now().time_since_epoch())).count();
}

long long cur_time_micros(){
  return (std::chrono::duration_cast< microseconds >(high_resolution_clock::now().time_since_epoch())).count();
}
string graph_extension = ".dbg";
string contig_extension = ".fasta";

struct parameters_t {
  std::string input_filename = "";
  std::string output_prefix = "";
};

void parse_arguments(int argc, char **argv, parameters_t & params);
void parse_arguments(int argc, char **argv, parameters_t & params)
{
  TCLAP::CmdLine cmd("Cosmo Copyright (c) Alex Bowe (alexbowe.com) 2014", ' ', VERSION);
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input",
            ".dbg file (output from cosmo-build).", true, "", "input_file", cmd);
  string output_short_form = "output_prefix";
  TCLAP::ValueArg<std::string> output_prefix_arg("o", "output_prefix",
            "Output prefix. Contigs will be written to [" + output_short_form + "]" + contig_extension + ". " +
            "Default prefix: basename(input_file).", false, "", output_short_form, cmd);
  cmd.parse( argc, argv );

  // -d flag for decompression to original kmer biz
  params.input_filename  = input_filename_arg.getValue();
  params.output_prefix   = output_prefix_arg.getValue();
}

int main(int argc, char* argv[]) {

  string index_filename = argv[1];
  string query_filename = argv[2];

  /*parameters_t p;
  parse_arguments(argc, argv, p);

  auto base_path = boost::filesystem::path(p.input_filename).parent_path().string();
  auto base_name = base_path + "/" + boost::filesystem::path(p.input_filename).stem().string();
  COSMO_LOG(debug) << base_name;
*/
  // TO LOAD:
  debruijn_graph_shifted<> g;
  load_from_file(g, index_filename);

  auto dbg_size = size_in_mega_bytes(g);
  cerr << "k             : " << g.k << endl;
  cerr << "num_nodes()   : " << g.num_nodes() << endl;
  cerr << "num_edges()   : " << g.num_edges() << endl;
  cerr << "W size        : " << size_in_mega_bytes(g.m_edges) << " MB" << endl;
  cerr << "L size        : " << size_in_mega_bytes(g.m_node_flags) << " MB" << endl;
  cerr << "DBG size      : " << dbg_size << " MB" << endl;
  cerr << "Bits per edge : " << bits_per_element(g) << " Bits" << endl;

  SeqIO::Reader<> in(query_filename);
  vector<string> queries;
  std::string str;
  while(true){
    long long len = in.get_next_read_to_buffer();
    if(len == 0) break;
    queries.push_back(std::string(in.read_buf));
  }

  cerr << "running " << queries.size() << " queries" << endl;
  // auto begin = chrono::high_resolution_clock::now();
  long long total_micros = 0;
  for (string s : queries) {
    long long t0 = cur_time_micros();
    auto res = g.index(s.begin());
    auto result = *res;
    total_micros += cur_time_micros() - t0;
    // if (s.substr(0, s.size() - 1) != g.node_label(g._edge_to_node(get<0>(result))))
    //   cerr << s.substr(0, s.size() - 1) << " " << g.node_label(g._edge_to_node(get<0>(result))) << endl;
  }
  
  cout << "Total query time us/kmer without I/O: " << (double)total_micros / queries.size() << endl;
  return 0;
}

