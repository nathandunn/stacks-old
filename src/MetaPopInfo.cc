#include <fstream>
#include <dirent.h>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <vector>

#include "constants.h"
#include "input.h"
#include "MetaPopInfo.h"

using std::ifstream;
using std::cerr;
using std::exception;

const string MetaPopInfo::Pop::default_name = "defaultpop";
const string MetaPopInfo::Group::default_name = "defaultgrp";

void MetaPopInfo::reset_sample_map() {
    sample_indexes_.clear();
    for (size_t i = 0; i < samples_.size(); ++i)
        sample_indexes_.insert( {samples_[i].name, i} );
}

void MetaPopInfo::reset_pop_map() {
    pop_indexes_.clear();
    for (size_t i = 0; i < pops_.size(); ++i)
        pop_indexes_.insert( {pops_[i].name, i} );
}

void MetaPopInfo::reset_group_map() {
    group_indexes_.clear();
    for (size_t i = 0; i < groups_.size(); ++i)
        group_indexes_.insert( {groups_[i].name, i} );
}

bool MetaPopInfo::init_popmap(const string& pmap_path, const string& dir_path) {

    ifstream fh(pmap_path.c_str(), ifstream::in);
    if (fh.fail())
        return false;

    size_t p = 0; // pop index counter
    size_t g = 0; // group index counter

    char line[max_len];
    memset(line, '\0', max_len);
    vector<string> parts;
    while (fh.getline(line, max_len)) {
        size_t len = strlen(line);

        // Skip empty lines and comments
        if (len == 0 || line[0] == '#')
            continue;

        // Check for Windows line endings
        if (line[len - 1] == '\r') {
            line[len - 1] = '\0';
            len -= 1;
        }

        // Parse the contents, we expect:
        // <file name> tab <population string> [tab <group string>]
        parse_tsv(line, parts);

        if (parts.size() < 2 || parts.size() > 3) {
            cerr << "Error: Malformed population map. File '" << pmap_path << "', line :\n" << line << "\n";
            throw exception();
        }

        //
        // Process the sample name.
        //

        samples_.push_back(Sample(parts[0]));

        // If the [dir_path] argument was given :
        // Check that "DIR/SAMPLENAME.matches.tsv(.gz)" exists.
        if(not dir_path.empty()) {
            bool exists = false;

            string filename = dir_path + samples_.back().name + ".matches.tsv";
            ifstream test_fh (filename);
            if (test_fh.good()) {
                exists = true;
            } else {
#ifdef HAVE_LIBZ
                gzFile gz_test_fh = gzopen((filename + "gz").c_str(), "rb");
                if (gz_test_fh) {
                    exists = true;
                }
                gzclose(gz_test_fh);
#endif
            }

            if (!exists) {
                cerr << " Unable to find '" << filename << "(.gz)', excluding this sample from the analysis.\n";
                samples_.pop_back();
                continue;
            }
        }

        //
        // Process the population field.
        //

        pair<map<string,size_t>::iterator, bool> pop_ins = pop_indexes_.insert( {parts[1], p} );
        size_t pop_index = pop_ins.first->second;

        samples_.back().pop = pop_index; // Set the sample's population index.
        if (pop_ins.second) {
            // Unknown pop
            pops_.push_back(Pop(parts[1]));
            ++p;
        }

        //
        // Process the group field, if any
        //

        if (parts.size() == 3) {
            pair<map<string,size_t>::iterator, bool> grp_ins = group_indexes_.insert( {parts[2], g} );
            size_t grp_index = grp_ins.first->second;

            if (grp_ins.second) {
                // Unknown group
                groups_.push_back(Group(parts[2]));
                ++g;
            }

            // If the current pop did not have a group yet, indicate
            // it. If it had one, check that it is the same.
            if (pops_[pop_index].group == size_t(-1)) {
                pops_[pop_index].group = grp_index;
                groups_[grp_index].pops.push_back(pop_index);
            } else if (pops_[pop_index].group != grp_index) {
                cerr << "Warning: In population map file '"
                     << pmap_path << "': population '"
                     << pops_[pop_index].name << "' belongs to two groups, '"
                     << groups_[pops_[pop_index].group].name << "' and '"
                     << groups_[grp_index].name << "'. Ignoring the latter one.\n";
            }
        }
    }
    if (samples_.size() == 0)
        return false;

    //
    // Check that all the populations are in a group. Put
    // populations that do not have a group in a default
    // one.
    //

    bool missing_group = false;
    if (not groups_.empty()) {
        for (vector<Pop>::iterator p = pops_.begin(); p != pops_.end(); ++p) {
            if (p->group == size_t(-1)) {
                cerr << "Warning: Population '" << p->name
                << "' did not have a group, adding it to '"
                << Group::default_name << "'.\n";
                missing_group = true;
            }
        }
    } else {
        missing_group = true;
    }
    if (missing_group) {
        groups_.push_back(Group(Group::default_name));
        g = groups_.size()-1;
        group_indexes_.insert( {Group::default_name, g} );
        for (size_t p = 0; p < pops_.size(); ++p) {
            if (pops_[p].group == size_t(-1)) {
                pops_[p].group = g;
                groups_[g].pops.push_back(p);
            }
        }
    }

    //
    // Sort the samples. Determine the first/last indexes for each population.
    //

    sort(samples_.begin(), samples_.end());
    reset_sample_map();

    size_t curr_pop = 0;
    pops_[curr_pop].first_sample = 0;
    for (size_t s = 1; s < samples_.size(); ++s) {
        if (samples_[s].pop != curr_pop) {
            pops_[curr_pop].last_sample = s-1;
            ++curr_pop;
            pops_[curr_pop].first_sample = s;
        }
    }
    pops_[curr_pop].last_sample = samples_.size()-1;

    return true;
}

bool MetaPopInfo::init_directory(const string& dir_path) {

    //
    // Find all sample names.
    //

    DIR* dir = opendir(dir_path.c_str());
    if (dir == NULL) {
        cerr << "Unable to open directory '" << dir_path << "' for reading.\n";
        throw exception();
    }
    dirent *direntry;
    while ((direntry = readdir(dir)) != NULL) {
        string filename = direntry->d_name;

        if (filename == "." || filename == ".." || filename.substr(0, 6) == "batch_")
            continue;

        size_t pos = filename.rfind(".tags.tsv");
        if (pos == string::npos)
            pos = filename.rfind(".tags.tsv.gz");

        if (pos != string::npos)
            samples_.push_back(Sample(filename.substr(0, pos)));
    }
    closedir(dir);

    if (samples_.size() == 0)
        return false;

    sort(samples_.begin(), samples_.end());

    //
    // Create a default group and a default population.
    //

    groups_.push_back(Group(Group::default_name));
    groups_.back().pops.push_back(0);

    pops_.push_back(Pop(Pop::default_name));
    pops_.back().first_sample = 0;
    pops_.back().last_sample = samples_.size()-1;
    pops_.back().group = 0;

    for (vector<Sample>::iterator s = samples_.begin(); s != samples_.end(); ++s)
        s->pop = 0;

    // Set the support members.
    reset_sample_map();
    reset_pop_map();
    reset_group_map();

    return true;
}

void MetaPopInfo::purge_samples(const vector<size_t>& rm_samples) {

    // Remove these samples from [samples_].
    for (vector<size_t>::const_iterator s = rm_samples.begin(); s != rm_samples.end(); ++s) {
        samples_.at(*s).name.clear(); // Mark the sample for removal.
    }
    samples_.erase(
            remove_if(samples_.begin(), samples_.end(),
                    [](Sample& s) {return s.name.empty();} ),
            samples_.end());

    // Update the indexes of the populations.
    for (vector<Pop>::iterator p = pops_.begin(); p != pops_.end(); ++p) {
        for (vector<size_t>::const_reverse_iterator rm_sample = rm_samples.rbegin(); rm_sample != rm_samples.rend(); ++rm_sample) {
            if (p->first_sample > *rm_sample) // n.b. ">"
                --p->first_sample;
            if (p->last_sample >= *rm_sample) // n.b. ">=". Thus if the population becomes
                                              // empty, [first_sample] will be past [last_sample].
                --p->last_sample;
        }
    }

    auto pop_is_empty = [](Pop& p){return (p.first_sample > p.last_sample);};

    // Remove the empty populations from [groups_].
    for(vector<Group>::iterator group = groups_.begin(); group != groups_.end(); ++group)
        group->pops.erase(
                remove_if(group->pops.begin(), group->pops.end(),
                        [this,&pop_is_empty](size_t p){return pop_is_empty(pops_[p]);}),
                group->pops.end());

    // Remove the empty populations from [pops_].
    pops_.erase(
            remove_if(pops_.begin(), pops_.end(), pop_is_empty),
            pops_.end());

    // Remove empty groups from [groups_].
    groups_.erase(
            remove_if(groups_.begin(), groups_.end(),
                    [](Group& g) {return g.pops.empty();}),
            groups_.end());

    // Update the support members.
    reset_sample_map();
    reset_pop_map();
    reset_group_map();
    reset_sample_id_map();
}

void MetaPopInfo::reset_sample_id_map() {
    sample_indexes_by_id_.clear();
    for (size_t i = 0; i < samples_.size(); ++i)
        sample_indexes_by_id_[samples_[i].id] = i;
}

void MetaPopInfo::fill_files(vector<pair<int, string> >& files) const {
    files.clear();
    for (vector<Sample>::const_iterator sample = samples_.begin(); sample != samples_.end(); ++sample)
        files.push_back( {sample->pop+1, sample->name} ); // i+1
}

void MetaPopInfo::fill_sample_ids(vector<int>& sample_ids) const {
    sample_ids.clear();
    for (vector<Sample>::const_iterator sample = samples_.begin(); sample != samples_.end(); ++sample)
        sample_ids.push_back(sample->id);
}

void MetaPopInfo::fill_samples(map<int, string>& samples) const {
    samples.clear();
    for (vector<Sample>::const_iterator sample = samples_.begin(); sample != samples_.end(); ++sample)
        samples.insert({sample->id, sample->name});
}

void MetaPopInfo::fill_pop_key(map<int, string>& pop_key) const {
    pop_key.clear();
    for (size_t i = 0; i < pops_.size(); ++i)
        pop_key.insert( {i+1, pops_[i].name} ); // i+1 (ids are indexes shifted by 1)
}

void MetaPopInfo::fill_pop_indexes(map<int, pair<int, int> >& pop_indexes) const {
    pop_indexes.clear();
    for (size_t i = 0; i < pops_.size(); ++i)
        pop_indexes.insert( {i+1, {pops_[i].first_sample, pops_[i].last_sample}} ); // i+1 (ids are indexes shifted by 1)
}

void MetaPopInfo::fill_grp_key(map<int, string>& grp_key) const {
    grp_key.clear();
    for (size_t i = 0; i < groups_.size(); ++i)
        grp_key.insert({i+1, groups_[i].name}); // i+1 (ids are indexes shifted by 1)
}

void MetaPopInfo::fill_grp_members(map<int, vector<int> >& grp_members) const {
    grp_members.clear();
    for (size_t i = 0; i < groups_.size(); ++i) {
        vector<int>& pop_ids = grp_members.insert( {i+1, vector<int>()} ).first->second; // i+1 (ids are indexes shifted by 1)
        for(vector<size_t>::const_iterator p = groups_[i].pops.begin(); p != groups_[i].pops.end(); ++p)
            pop_ids.push_back(*p+1); // p+1 (ids are indexes shifted by 1)
    }
}
