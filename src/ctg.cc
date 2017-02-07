#include <string>

#include <gatb/gatb_core.hpp>

using namespace std;

int main(int argc, char* argv [])
{
    // TODO Initialize graph?

    // Read the input sequences.
    {
        const string reads_path = "./reads.fq.gz";
        IBank * inbank = Bank::open(reads_path);
        Iterator<Sequence>* it = inbank->iterator();
        for (it->first(); !it->isDone(); it->next()) {
            Sequence& seq = it->item();
            // Print the sequence.
            cout << " [ " << seq.getDataSize() << " ] " << seq.getComment() << endl;
            cout << seq.toString() << endl ;

            // TODO Build graph?
        }
        delete it;
        delete inbank;
    }

    // TODO Process graph?
}

/*
 * // More code examples -- c.f. Tutorial "GATB Programming Day"
 *
 * // Writing sequence files.
 * IBank* outBank = new BankFasta("path") ;
 * for (...) {Sequence& seq = ...; outBank->insert(seq);}
 * outBank->flush();
 *
 * // Manipulating Kmers.
 * // (Also `KmerCanonical`.)
 * Kmer<span>::ModelDirect model (kmerSize);
 * Kmer<span>::ModelCanonical modelc (kmerSize);
 * Kmer <span>::ModelDirect::Iterator itKmer (model);
 * itKmer->value() // forward(), revcomp()
 *
 * // Graph creation: Graph::create or Graph::load, e.g.
 * Graph graph = Graph::create("-in %s -abundance-min %d -kmer-size %d -debloom original", reads_path, 3, kmerSize);
 *
 * // Iterating over nodes -- c.f. "Overview"
 *
 * Graph::Iterator<Node> it = graph.iterator<Node>();
 * for (it.first(); !it.isDone(); it.next()) {
 *     std::cout << graph.toString (it.item()) << std::endl;
 * }
 */
