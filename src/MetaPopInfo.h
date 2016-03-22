/*
 * todo Ask @Julian, are you OK with the [variable] notation for comments ?
 */
#ifndef METAPOPINFO_H
#define METAPOPINFO_H

#include <string>
#include <vector>
#include <map>

using std::size_t;
using std::pair;
using std::vector;
using std::string;
using std::map;

/*
 * MetaPopInfo
 * Class for reprensenting a metapopulation : its individuals/samples,
 * populations, and groups of populations.
 */
class MetaPopInfo {
public:
    typedef size_t sample_index; // For indexes in [samples_].
    typedef size_t pop_index;
    typedef size_t group_index;

    struct Sample {
        pop_index pop;
        string prefix;
        size_t id; // Sample ids, as present in the matches files.
    };
    struct Pop {
        group_index group;
        sample_index first_sample;
        sample_index last_sample;
        string name;
    };
    struct Group {
        vector<pop_index> pops;
        string name;
    };

private:
    vector<Sample> samples_; //n.b. Samples must be grouped by population.
    vector<Pop> pops_;
    vector<Group> groups_;

    map<size_t,sample_index> sample_indexes_by_id_; // Links an id with an index in [samples_].
    map<string,sample_index> sample_indexes_by_name_; // Links a name with an index in [samples_].
    map<string,pop_index> pop_indexes_; // same, for populations
    map<string,group_index> group_indexes_; // same, for groups

public:
    /*
     * Generate a description of the population according to the input :
     * -- If a popmap file was provided, load its contents.
     * -- Otherwise, browse the directory.
     */
    void init(const string& dir_path, const string& popmap_relpath = string()); //todo implement it.

    // Add an sstacks id to a sample
    void set_sample_id(const sample_index i, const size_t id) {samples_.at(i).id = id;}

    // Access to the information
    const vector<Sample>& samples() const {return samples_;}
    const vector<Pop>& pops() const {return pops_;}
    const vector<Group>& groups() const {return groups_;}

    // Obtain the indexes corresponding to a particular id or name.
    size_t get_sample_index(const string& name) const {return sample_indexes_by_name_.at(name);}
    size_t get_sample_index(const size_t& id) const {return sample_indexes_by_id_.at(id);}
    size_t get_pop_index(const string& name) const {return pop_indexes_.at(name);}
    size_t get_group_index(const string& name) const {return group_indexes_.at(name);}

    namespace backcompat {
    void fill_files(vector<pair<int, string> >&) const; // vector of (pop_index, prefix) pairs
    void fill_samples(map<int, string>&) const; // map of sample_id : sample_name
    void fill_sample_ids(vector<int>&) const; // vector of sample_id's.
    void fill_pop_key(map<int, string>&) const; // map of pop_index : pop_name
                                                // n.b. pop_indexes start at 1
    void fill_pop_indexes(map<int, pair<int, int> >&) const; // map of pop_index : (sample_index, sample_index)
    void fill_grp_key(map<int, string>&) const; // map of group_index : group_name
                                                // n.b. group indexes start at 1
    void fill_grp_members(map<int, vector<int> >) const; // map of group_index : pop_indexes
    }
};

#endif // METAPOPINFO_H
