/*
 * todo Ask @Julian, are you OK with the [variable] notation for comments ?
 */
#ifndef METAPOPULATION_DESCRIPTION_H
#define METAPOPULATION_DESCRIPTION_H

#include <string>
#include <vector>
#include <map>
using namespace std;

class Metapopulation_description;
typedef Metapopulation_description Mpopdesc;

class Metapopulation_description {
public:
    typedef size_t sample_index; // For indexes in [samples_].
    typedef size_t pop_index;
    typedef size_t group_index;

    struct Sample {
        pop_index pop;
        string name;
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

    map<string,sample_index> sample_indexes_; // Links a name with an index in [samples_].
    map<string,pop_index> pop_indexes_;
    map<string,group_index> group_indexes_;

public:
    void load(const string& path); // Loads a popmap file

    const vector<Sample>& samples() const {return samples_;}
    const vector<Pop>& pops() const {return pops_;}
    const vector<Group>& groups() const {return groups_;}

    size_t get_sample_index(const string& name) const {return sample_indexes_.at(name);}
    size_t get_pop_index(const string& name) const {return pop_indexes_.at(name);}
    size_t get_group_index(const string& name) const {return group_indexes_.at(name);}

};

void Metapopulation_description::load(const string& path) {} // todo

#endif // METAPOPULATION_DESCRIPTION_H
