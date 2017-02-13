#include <vector>
#include <string>

#define WITH_LAMBDA_EXPRESSION
#include <gatb/gatb_core.hpp>

#include <gatb/debruijn/impl/Simplifications.hpp>
#include <gatb/debruijn/impl/GraphUnitigs.hpp>

using namespace std;

const string reads_path = "reads.fq.gz";
const string out_path = "out.fa";
const size_t min_cov = 2;
const size_t klen = 20;
const bool remove_tips = false;
const bool remove_bulges = false;
const bool remove_ec = false;

struct Contig {
    NodeGU starting_node;
    string seq;
    float cov;
    bool isolated_left;
    bool isolated_right;
};

//
// Classes to use with Integer::apply<Functor, Params>(size_t klen, Params params);
//

struct AssemblyArgs {
    const vector<string>& reads;
    size_t klen;
    size_t min_cov;
    vector<Contig>& contigs;
};

template<size_t span>
struct AssemblyFunctor {
    void operator() (AssemblyArgs& p);
};

//
// main
//

int main(int argc, char* argv [])
{
    // Read the sequences into a vector<string>
    vector<string> reads;
    {
        IBank* ibank = Bank::open (reads_path);
        Iterator<Sequence>* s = ibank->iterator();
        for(s->first(); !s->isDone(); s->next())
            reads.push_back(s->item().toString());
    }

    // Do the assembly
    vector<Contig> contigs;
    AssemblyArgs args = {reads, klen, min_cov, contigs};
    Integer::apply<AssemblyFunctor, AssemblyArgs>(args.klen, args);

    // Write the resulting contigs
    cout << "args.contigs->size(): " << args.contigs.size() << endl;
    ofstream ofile (out_path);
    ofile << fixed << setprecision(2);
    for(size_t i=0; i<args.contigs.size(); ++i) {
        Contig& c = args.contigs.at(i);
        ofile << ">NODE_"<< i+1 << "_length_" << c.seq.size() << "_cov_" << c.cov << "_ID_" << i << "\n"
              << c.seq << "\n";
    }
}

template<size_t span>
void AssemblyFunctor<span>::operator() (AssemblyArgs& args) {
    typedef GraphUnitigsTemplate<span> Graph;

    // Initialize the graph
    BankStrings* bank = new BankStrings(args.reads);
    Graph graph = Graph::create(bank, "-kmer-size %d -abundance-min %d", args.klen, args.min_cov);
    graph.precomputeAdjacency();

    // Perform simplifications
    Simplifications<Graph,NodeGU,EdgeGU> simplifications(graph, 1, false); // one core, quiet
    simplifications._doTipRemoval = remove_tips;
    simplifications._doBulgeRemoval = remove_bulges;
    simplifications._doECRemoval = remove_ec;
    simplifications.simplify();

    // Retrieve the resulting contigs
    GraphIterator<NodeGU> node_it = graph.iterator();
    for(node_it.first(); !node_it.isDone(); node_it.next()) {
        NodeGU node = node_it.item();
        if (graph.unitigIsMarked(node) || graph.isNodeDeleted(node))
            continue;

        args.contigs.push_back(Contig());
        Contig& c = args.contigs.back();
        c.starting_node = node;
        c.seq = graph.simplePathBothDirections(c.starting_node, c.isolated_left, c.isolated_right, true, c.cov);
    }
}
