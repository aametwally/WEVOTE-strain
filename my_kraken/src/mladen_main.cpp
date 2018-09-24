#include "kraken_headers.hpp"
#include "krakendb.hpp"
#include "krakenutil.hpp"
#include "krakendb.hpp"
#include "quickfile.hpp"



using namespace std;
using namespace kraken;

//Current goal:
//1. open fasta file
//2. parse fasta contents and loop through file
//3. take dnasequence, and run through kmer stuff (get this from kraken)
//


const size_t DEF_WORK_UNIT_SIZE = 500000;

int Num_threads = 1;
string Nodes_filename;
string DB_filename="/export/home/mrasic2/ecoliDB/database.kdb";
string Index_filename="/export/home/mrasic2/ecoliDB/database.idx";
bool Quick_mode = false;
bool Fastq_input = false;
bool Fastq_output = false;
bool Paired_input = false;
bool Print_classified = false;
bool Print_unclassified = false;
bool Print_kraken = true;
bool Populate_memory = false;
bool Only_classified_kraken_output = false;
uint32_t Minimum_hit_count = 1;
map<uint32_t, uint32_t> Parent_map;
KrakenDB Database;
string Classified_output_file, Unclassified_output_file, Kraken_output_file;
string Output_format;
ostream *Classified_output;
ostream *Classified_output2;
ostream *Unclassified_output;
ostream *Unclassified_output2;
ostream *Kraken_output;
size_t Work_unit_size = DEF_WORK_UNIT_SIZE;

uint64_t total_classified = 0;
uint64_t total_sequences = 0;
uint64_t total_bases = 0;


void classify_sequence(string &dna, ostringstream &koss);
string hitlist_string(vector<uint32_t> &taxa, vector<uint8_t> &ambig);   

int main(int argc, char **argv){
    string current_line;
    string seq;
    string dnaseq;
    vector<string> dna_list;
    ostringstream k1,c1,u1,c2,u2;

    ifstream readsFile( argv[1] ); //read the file in argv1
    if (readsFile.is_open())
    {

	//OPEN file and store each DNAseq into the vector<string> called dna_list
        while ( getline (readsFile,current_line) )
        {
		if(current_line.size() > 0 && current_line[0] != '>')
		{
			dna_list.push_back(current_line);
		}
        }

    }


    QuickFile db_file;
    db_file.open_file(DB_filename);

    Database = KrakenDB(db_file.ptr());
    KmerScanner::set_k(Database.get_k());

    
    QuickFile idx_file;
    idx_file.open_file(Index_filename);
    KrakenDBIndex db_index(idx_file.ptr());
    Database.set_index(&db_index);


    if (! Kraken_output_file.empty()) {                          
        if (Kraken_output_file == "-")
            Print_kraken = false;
    else
        Kraken_output = new ofstream(Kraken_output_file.c_str());
    }
    else

        Kraken_output = &cout;


    k1.str("");
    c1.str("");
    u1.str("");
    c2.str("");
    u2.str("");



    for (size_t j=0; j < 300; j++)
        classify_sequence(dna_list[j],k1);

        (*Kraken_output) << k1.str();
    total_sequences += dna_list.size();

        if (isatty(fileno(stderr)))
            cerr << "hi";

  if (isatty(fileno(stderr)))
    cerr << "\r";

  fprintf(stderr, "  %llu sequences classified (%.2f%%)\n",
          (unsigned long long) total_classified, total_classified * 100.0 / total_sequences);
  fprintf(stderr, "  %llu sequences unclassified (%.2f%%)\n",
          (unsigned long long) (total_sequences - total_classified),
          (total_sequences - total_classified) * 100.0 / total_sequences);

    return 0;
}








void classify_sequence(string &dna, ostringstream &koss)
{
    vector<uint32_t> taxa;
    vector<uint8_t> ambig_list;
    map<uint32_t, uint32_t> hit_counts;
    uint64_t *kmer_ptr;
    uint32_t taxon = 0;
    uint32_t hits = 0;  // only maintained if in quick mode
    uint64_t current_bin_key;
    int64_t current_min_pos = 1;
    int64_t current_max_pos = 0;

  if (dna.size() >= Database.get_k()) {
    KmerScanner scanner(dna);
    while ((kmer_ptr = scanner.next_kmer()) != NULL) {
      taxon = 0;
      if (scanner.ambig_kmer()) {
        ambig_list.push_back(1);
      }
      else {
        ambig_list.push_back(0);
        uint32_t *val_ptr = Database.kmer_query(
                              Database.canonical_representation(*kmer_ptr),
                              &current_bin_key,
                              &current_min_pos, &current_max_pos
                            );
        taxon = val_ptr ? *val_ptr : 0;
        if (taxon) {
          hit_counts[taxon]++;
          if (Quick_mode && ++hits >= Minimum_hit_count)
            break;
        }
      }
      taxa.push_back(taxon);
    }
  }


  uint32_t call = 0;
  call = resolve_tree(hit_counts, Parent_map);

  if (call)
    total_classified++;

  if (! Print_kraken)
    return;

  if (call) {
    koss << "C\t";
  }
  else {
    koss << "U\t";
  }
  //koss << dna.id << "\t" << call << "\t" << dna.size() << "\t";
  //replacing it with my version (no dna.id saved)
  koss << call << "\t" << dna.size() << "\t";

  if (taxa.empty())
      koss <<"0:0";
  else
      koss << hitlist_string(taxa,ambig_list);

  koss << endl;

}

string hitlist_string(vector<uint32_t> &taxa, vector<uint8_t> &ambig)                                              
{
  int64_t last_code;
  int code_count = 1;
  ostringstream hitlist;

  if (ambig[0])   { last_code = -1; }
  else            { last_code = taxa[0]; }

  for (size_t i = 1; i < taxa.size(); i++) {
    int64_t code;
    if (ambig[i]) { code = -1; }
    else          { code = taxa[i]; }

    if (code == last_code) {
      code_count++;
    }
    else {
      if (last_code >= 0) {
        hitlist << last_code << ":" << code_count << " ";
      }
      else {
        hitlist << "A:" << code_count << " ";
      }
      code_count = 1;
      last_code = code;
    }
  }
  if (last_code >= 0) {
    hitlist << last_code << ":" << code_count;
  }
  else {
    hitlist << "A:" << code_count;
  }
  return hitlist.str();
}










