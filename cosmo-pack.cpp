#include <iterator>
#include <typeinfo>
#include <fstream>
#include <libgen.h>
#include <stxxl.h>

#include <stxxl/bits/containers/sorter.h>
#include <stxxl/bits/parallel.h>
#include <boost/range/adaptors.hpp>
#include <boost/range/adaptor/uniqued.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/function_input_iterator.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/core.hpp>
#include <boost/log/expressions.hpp>

// TCLAP
#include "tclap/CmdLine.h"

#include "utility.hpp"
#include "kmer.hpp"
#include "io.hpp"
#include "sort.hpp"
#include "dummies.hpp"

//#define BOOST_LOG_DYN_LINK 1
#define COSMO_LOG BOOST_LOG_TRIVIAL
#define STXXL_FORCE_VERBOSE_LEVEL -10

#if K_LEN <= 32
typedef uint64_t kmer_t;
#elif K_LEN <= 64
typedef uint128_t kmer_t;
#endif // add some static error handling to this
typedef std::pair<kmer_t, size_t>  dummy_t;
const static string extension = ".packed";

const size_t block_size = 2 * 1024 * 1024; // 512 MB -> safe value for cuda?

typedef node_less<kmer_t> node_comparator_type;
typedef kmer_less<kmer_t> edge_comparator_type;
typedef colex_dummy_less<dummy_t> dummy_comparator_type;

// TODO: Work out how to hijack stxxl inmemory sort
/*
namespace stxxl {
namespace potentially_parallel {
template <class Iterator>
void sort(Iterator a, Iterator b, node_comparator_type cmp) {
  cerr << "Foo" << endl;
  assert(0);
}
}
}
*/
typedef stxxl::sorter<kmer_t, node_comparator_type, block_size> node_sorter_type;
typedef stxxl::sorter<kmer_t, edge_comparator_type, block_size> edge_sorter_type;
typedef stxxl::sorter<dummy_t, dummy_comparator_type, block_size> dummy_sorter_type;

// TODO: output dummy positions instead
typedef struct p
{
  //bool ascii = false;
  std::string input_filename = "";
  std::string output_prefix = "";
  size_t k = 0;
  size_t m = 0;
} parameters_t;

void parse_arguments(int argc, char **argv, parameters_t & params);
void parse_arguments(int argc, char **argv, parameters_t & params)
{
  TCLAP::CmdLine cmd("Cosmo Copyright (c) Alex Bowe (alexbowe.com) 2014", ' ', VERSION);
  /* // Add this option after refactoring the visitors (for now just compile with DEBUG if you want printed edges)
  TCLAP::SwitchArg ascii_arg("a", "ascii",
            "Outputs *full* edges (instead of just last nucleotide) as ASCII.",
            cmd, false);
  */
  TCLAP::UnlabeledValueArg<std::string> input_filename_arg("input",
            "Input file. Currently only supports DSK's binary format (for k<=64).", true, "", "input_file", cmd);
  string output_short_form = "output_prefix";
  TCLAP::ValueArg<size_t> kmer_length_arg("k", "kmer_length", "Length of edges (node is k-1).", true, 0, "length", cmd);
  size_t default_mem = 4 * 1024;
  TCLAP::ValueArg<size_t> mem_size_arg("m", "mem_size", "Internal memory to use (MB).", false, default_mem, "length", cmd);

  TCLAP::ValueArg<std::string> output_prefix_arg("o", "output_prefix",
            "Output prefix. Results will be written to [" + output_short_form + "]" + extension + ". " +
            "Default prefix: basename(input_file).", false, "", output_short_form, cmd);
  cmd.parse( argc, argv );
  //params.ascii         = ascii_arg.getValue();
  params.input_filename  = input_filename_arg.getValue();
  params.k               = kmer_length_arg.getValue();
  params.m               = mem_size_arg.getValue() * MB_TO_BYTES;
  params.output_prefix   = output_prefix_arg.getValue();
}

int main(int argc, char* argv[])
{
  using namespace boost::adaptors;

  // Parameter extraction
  parameters_t params;
  parse_arguments(argc, argv, params);
  stxxl::internal_size_type M = params.m;
  std::string file_name = params.input_filename;
  char * base_name = basename(const_cast<char*>(file_name.c_str()));
  size_t k = params.k;

  // TODO: Logging/verbosity flags
  // TODO: Config file
  // TODO: ASCII outputter (first and last symbols, option for whole kmer)
  // TODO: stdin input
  boost::log::core::get()->set_filter (
    boost::log::trivial::severity != boost::log::trivial::debug
  );
  /*
  COSMO_LOG(trace)   << "A trace severity message";
  COSMO_LOG(debug)   << "A debug severity message";
  COSMO_LOG(info)    << "An informational severity message";
  COSMO_LOG(warning) << "A warning severity message";
  COSMO_LOG(error)   << "An error severity message";
  COSMO_LOG(fatal)   << "A fatal severity message";
  */


  // Check k value supported
  if (k > K_LEN) {
    COSMO_LOG(error) << "This version only supports k <= " << K_LEN << ". Try recompiling.";
    exit(1);
  }

  // TODO: Test sort input size (with sorter ctor, overload sort func)
  // block size vs creator size?

  // Open the file
  stxxl::syscall_file in_file(file_name, stxxl::file::DIRECT | stxxl::file::RDONLY);
  //stxxl::syscall_file out_file(file_name + ".boss",
  //stxxl::file::DIRECT | stxxl::file::RDWR | stxxl::file::CREAT);
  typedef stxxl::vector<kmer_t, 1, stxxl::lru_pager<8>, block_size> vector_type;
  vector_type in_vec(&in_file); 
  //vector_type kmers(&out_file);
  vector_type kmers;
  kmers.resize(in_vec.size()*2); // x2 for revcomps
  vector_type kmers_b;
  kmers_b.resize(in_vec.size()*2); // x2 for revcomps
  node_sorter_type node_sorter(node_comparator_type(), M/2);
  edge_sorter_type edge_sorter(edge_comparator_type(), M/2);

  // Make conversion functors
  auto swap  = swap_gt<kmer_t>();
  auto revnt = reverse_nt<kmer_t>();
  auto rc    = reverse_complement<kmer_t>(k);

  // Convert to our format: reverse for colex ordering, swap g/t encoding (DSK)
  COSMO_LOG(trace) << "Creating runs...";
  // TODO: Parallelise? #pragma omp parallel for
  for (kmer_t record : vector_type::bufreader_type(in_vec)) {
    auto x = revnt(swap(record));
    auto y = rc(x);
    node_sorter.push(x);
    node_sorter.push(y);
    edge_sorter.push(x);
    edge_sorter.push(y);
  }
  COSMO_LOG(info) << "Added " << node_sorter.size()/2 << " kmers, and their reverse complements.";

  // TODO: replace internal sort with nvidia radix sort (cub/thrust)
  // TODO: test simulated B table vs second sort
  COSMO_LOG(trace) << "Merging runs...";
  node_sorter.sort();
  edge_sorter.sort();

  // TODO: make buffered reader around sorted stream instead
  // Or keep sorting two tables and stream to range for dummy edge finding (then rewind)
  
  COSMO_LOG(trace) << "Writing to temporary storage...";
  stxxl::stream::materialize(node_sorter, kmers.begin(), kmers.end());
  node_sorter.finish_clear();
  stxxl::stream::materialize(edge_sorter, kmers_b.begin(), kmers_b.end());
  edge_sorter.finish_clear();

  // TODO: read about STXXL_PARALLEL_MULTIWAY_MERGE and other defs
  // Find nodes that require incoming dummy edges
  /*
  vector_type::bufreader_type br(kmers);
  auto f= std::function<kmer_t(void)>([&]()->kmer_t{
    kmer_t temp = *br;
    ++br;
    return temp;
  });
  typedef boost::function_input_iterator<decltype(f), size_t> fiit;
  fiit foo_0 = fiit(f, (size_t)0);
  fiit foo_1 = fiit(f, kmers.size());
  auto foo = boost::make_iterator_range(foo_0, foo_1);
  */
  // TODO: Write sort stream (or buffered reader) to iterator adapter
  // TODO: Make iterator adapter for stream/bufreader (*, ++)
  
  
  stxxl::vector<dummy_t> incoming_dummies;
  #define CI_RANGE(x) (boost::make_iterator_range((x).cbegin(),(x).cend()))
  auto a = CI_RANGE(kmers);
  auto b = CI_RANGE(kmers_b);
  // TODO: function output to buffered writer or sorter, for dummy edges
  // or a b-tree if generating all shifts
  COSMO_LOG(trace) << "Searching for nodes requiring incoming dummy edges...";

  dummy_sorter_type dummy_sorter(dummy_comparator_type(), M);
  // TODO: redo this so we dont need the dummies (follow up with Travis)
  find_incoming_dummy_nodes<kmer_t>(a, b, k, [&](kmer_t x) {
    //dummy_sorter.push(dummy_t(x, k-1));
    for (int i=0; i<k-1; i++) {
       x <<= i*2;
       dummy_sorter.push(dummy_t(x,k-i-1));
    }
  });
  COSMO_LOG(info)  << "Added " << dummy_sorter.size() << " incoming dummy edges.";
  COSMO_LOG(trace) << "Sorting dummies...";
  dummy_sorter.sort();

  incoming_dummies.resize(dummy_sorter.size());
  stxxl::stream::materialize(dummy_sorter, incoming_dummies.begin(), incoming_dummies.end());

  // Make Outputter
  // TODO: Should probably do checking here when opening the file...
  string outfilename = (params.output_prefix == "")? base_name : params.output_prefix;
  ofstream ofs;
  ofs.open(outfilename + extension, ios::out | ios::binary);
  PackedEdgeOutputer out(ofs);
  #ifdef VAR_ORDER
  ofstream lcs;
  lcs.open(outfilename + extension + ".lcs", ios::out | ios::binary);
  #endif

  COSMO_LOG(trace) << "Merging dummies and outputting...";
  // Merge dummies and output
  size_t prev_k = 0;
  // TODO: rewrite to use bufreaders
  // TODO: **** ADD COLOUR PASSTHROUGH **** (check this in first)
  auto d = CI_RANGE(incoming_dummies);
  merge_dummies(a, b, d, k,
    [&](edge_tag tag, const kmer_t & x, size_t this_k, size_t lcs_len, bool first_end_node) {
      #ifdef VAR_ORDER
      out.write(tag, x, this_k, (lcs_len != k-1), first_end_node);
      char l(lcs_len);
      lcs.write((char*)&l, 1);
      #else
      out.write(tag, x, this_k, lcs_len, first_end_node);
      #endif
      prev_k = this_k;
      
      #ifdef VERBOSE // print each kmer to stderr for testing
      if (tag == out_dummy) cout << kmer_to_string(get_start_node(x), k-1, k-1) << "$";
      else                  cout << kmer_to_string(x, k, this_k);
      cout << " " << (lcs_len != k-1) << " " << lcs_len << " " << first_end_node << endl;
      #endif
    });

  
  out.close();
  #ifdef VAR_ORDER
  lcs.flush();
  lcs.close();
  #endif


  out.close();
  #ifdef VAR_ORDER
  lcs.flush();
  lcs.close();
  #endif

  // Write out counts and stuff
  uint64_t t_k(k); // make uint64_t just to make parsing easier
  // (can read them all at once and take the last 6 values)
  ofs.write((char*)&t_k, sizeof(uint64_t));
  ofs.flush();
  ofs.close();
  COSMO_LOG(trace) << "Done!";

  return 0;
}
