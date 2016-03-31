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
    struct Sample {
        string name;
        size_t pop;
        size_t id; // optional, deprecated

        Sample(const string& name) : name(name), pop(-1), id(-1) {}
        inline bool operator<(const Sample& other) const;
    };
    struct Pop {
        string name;
        size_t first_sample;
        size_t last_sample;
        size_t group;

        Pop(const string& name) : name(name), group(-1), first_sample(-1), last_sample(-1) {}
        static const string default_name;
    };
    struct Group {
        string name;
        vector<size_t> pops;

        Group(const string& name) : name(name) {}
        static const string default_name;
    };

private:
    vector<Sample> samples_; //n.b. Samples are sorted primarily by population index, and secondarily by name.
    vector<Pop> pops_;
    vector<Group> groups_;

    map<string,size_t> sample_indexes_; // Links a name with an index in [samples_].
    map<string,size_t> pop_indexes_;
    map<string,size_t> group_indexes_;
    void reset_sample_map(); // Resets [sample_indexes_].
    void reset_pop_map();
    void reset_group_map();

public:
    // Create the representation :
    // -- from a population map file. For consistency with the existing
    //    code, the existence of the "DIR/SAMPLENAME.matches.tsv(.gz)"
    //    files is checked if a [dir_path] argument is given.
    // -- or by browsing the directory for "*.tags.tsv(.gz)" files.
    bool init_popmap(const string& popmap_path, const string& dir_path = string());
    bool init_directory(const string& dir_path);

    // Removes samples from the metapopulation.
    // As samples, populations or groups are removed, the indexes of
    // the remaining ones change, but the order in which they appear
    // is preserved.
    void purge_samples(const vector<size_t>& samples);

    // Retrieve information.
    const vector<Sample>& samples() const {return samples_;}
    const vector<Pop>& pops() const {return pops_;}
    const vector<Group>& groups() const {return groups_;}

    size_t get_sample_index(const string& name) const {return sample_indexes_.at(name);}
    size_t get_pop_index(const string& name) const {return pop_indexes_.at(name);}
    size_t get_group_index(const string& name) const {return group_indexes_.at(name);}

private:
    map<size_t,size_t> sample_indexes_by_id_; // Links a sample id with an index in [samples_].

public:
    /*
     * Methods for backwards compatibility
     */

    // Sets the ID of a sample.
    void set_sample_id(size_t index, size_t id) {samples_.at(index).id = id; sample_indexes_by_id_[id] = index;}
    size_t get_sample_index(const size_t& id) const {return sample_indexes_by_id_.at(id);}

    // Resets the (sample_id : index) map. It is the caller's responsibility that the ids are unique.
    void reset_sample_id_map();

    void fill_files(vector<pair<int, string> >&) const;
    void fill_sample_ids(vector<int>&) const;
    void fill_samples(map<int, string>&) const;
    void fill_pop_key(map<int, string>&) const;
    void fill_pop_indexes(map<int, pair<int, int> >&) const;
    void fill_grp_key(map<int, string>&) const;
    void fill_grp_members(map<int, vector<int> >&) const;
};

bool MetaPopInfo::Sample::operator<(const Sample& other) const {
    if (pop == other.pop)
        return name < other.name;
    else
        return pop < other.pop;
}

#endif // METAPOPINFO_H
