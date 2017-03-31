#include "sql_utilities.h"

using namespace std;

void load_catalog_matches(string sample,  vector<CatMatch *> &matches) {
    CatMatch      *m;
    string         f;
    vector<string> parts;
    long int       line_num;
    ifstream       fh;
    gzFile         gz_fh;

    char *line      = (char *) malloc(sizeof(char) * max_len);
    int   size      = max_len;
    int   cnt       = 0;
    bool  gzip      = false;
    int   fh_status = 1;

    f = sample + ".matches.tsv";
    fh.open(f.c_str(), ifstream::in);
    if (fh.fail()) {
        //
        // Test for a gzipped file.
        //
        f = sample + ".matches.tsv.gz";
        gz_fh = gzopen(f.c_str(), "rb");
        if (!gz_fh)
            return;

        #if ZLIB_VERNUM >= 0x1240
        gzbuffer(gz_fh, libz_buffer_size);
        #endif
        gzip = true;
    }
    cerr << "  Parsing " << f.c_str() << "\n";

    line_num = 1;
    while (fh_status) {
        fh_status = (gzip == true) ? read_gzip_line(gz_fh, &line, &size) : read_line(fh, &line, &size);
        line_num++;

        if (!fh_status && strlen(line) == 0)
            continue;

        if (is_comment(line)) continue;

        parse_tsv(line, parts);

        cnt = parts.size();

        if (cnt != num_matches_fields && cnt != num_matches_fields - 1) {
            cerr << "Error parsing " << f.c_str() << " at line: " << line_num << ". (" << parts.size() << " fields).\n";
            throw exception();
        }

        m = new CatMatch;
        m->batch_id  = atoi(parts[1].c_str());
        m->cat_id    = atoi(parts[2].c_str());
        m->sample_id = atoi(parts[3].c_str());
        m->tag_id    = atoi(parts[4].c_str());
        m->haplotype = new char[parts[5].length() + 1];
        strcpy(m->haplotype, parts[5].c_str());
        m->depth     = atoi(parts[6].c_str());
        m->lnl       = is_double(parts[7].c_str());

        if (cnt == num_matches_fields && parts[8].length() > 0) {
            m->cigar = new char[parts[8].length() + 1];
            strcpy(m->cigar, parts[8].c_str());
        }

        matches.push_back(m);
    }

    if (gzip)
        gzclose(gz_fh);
    else
        fh.close();

    return;
}

int load_model_results(string sample,  map<int, ModRes *> &modres) {
    string         f;
    vector<string> parts;
    long int       line_num;
    ifstream       fh;
    gzFile   gz_fh;

    char *line      = (char *) malloc(sizeof(char) * max_len);
    int   size      = max_len;
    bool  gzip      = false;
    bool  open_fail  = false;
    int   fh_status = 1;

    //
    // Parse the models file (if it exists), otherwise parse the tag file to
    // pull in the model calls for each locus.
    //
    gzip      = false;
    fh_status = 1;
    line_num  = 1;

    f = sample + ".models.tsv";
    fh.open(f.c_str(), ifstream::in);
    if (fh.fail())
        open_fail = true;

    if (open_fail) {
        //
        // Test for a gzipped MODELs file.
        //
        f = sample + ".models.tsv.gz";
        gz_fh = gzopen(f.c_str(), "rb");
        if (!gz_fh) {
            open_fail = true;
        } else {
            open_fail = false;
            #if ZLIB_VERNUM >= 0x1240
            gzbuffer(gz_fh, libz_buffer_size);
            #endif
            gzip = true;
        }
    }

    if (open_fail) {
        //
        // Test for a TAGs file.
        //
        f = sample + ".tags.tsv";
        fh.open(f.c_str(), ifstream::in);
        if (fh.fail())
            open_fail = true;
        else
            open_fail = false;
    }

    if (open_fail) {
        //
        // Test for a gzipped TAGs file.
        //
        f = sample + ".tags.tsv.gz";
        gz_fh = gzopen(f.c_str(), "rb");
        if (!gz_fh) {
            open_fail = true;
        } else {
            open_fail = false;
            #if ZLIB_VERNUM >= 0x1240
            gzbuffer(gz_fh, libz_buffer_size);
            #endif
            gzip = true;
        }
    }

    if (open_fail) {
        cerr << " Unable to open '" << sample << "'\n";
        return 0;
    }

    cerr << "  Parsing " << f.c_str() << "\n";

    ModRes *mod;
    uint    tag_id, samp_id;

    while (fh_status) {
        fh_status = (gzip == true) ? read_gzip_line(gz_fh, &line, &size) : read_line(fh, &line, &size);
        line_num++;

        if (!fh_status && strlen(line) == 0)
            continue;
        if (is_comment(line)) continue;

        parse_tsv(line, parts);

        if (parts.size() != num_tags_fields) {
            cerr << "Error parsing " << f.c_str() << " at line: " << line_num << ". (" << parts.size() << " fields).\n";
            return 0;
        }

        //
        // Read the model sequence, a series of letters specifying if the model called a
        // homozygous base (O), a heterozygous base (E), or if the base type was unknown (U).
        //
        if (parts[6] != "model") continue;

        samp_id = atoi(parts[1].c_str());
        tag_id  = atoi(parts[2].c_str());
        mod     = new ModRes(samp_id, tag_id, parts[9].c_str());

        modres[tag_id] = mod;
    }

    if (gzip)
        gzclose(gz_fh);
    else
        fh.close();

    delete [] line;

    return 1;
}

int load_snp_calls(string sample,  map<int, SNPRes *> &snpres) {
    string         f;
    int id, samp_id;
    vector<string> parts;
    long int       line_num;
    SNP           *snp;
    SNPRes        *snpr;
    ifstream       fh;
    gzFile   gz_fh;

    char *line      = (char *) malloc(sizeof(char) * max_len);
    int   size      = max_len;
    bool  gzip      = false;
    int   fh_status = 1;

    //
    // Parse the SNP file
    //
    f = sample + ".snps.tsv";
    fh.open(f.c_str(), ifstream::in);
    if (fh.fail()) {
        //
        // Test for a gzipped file.
        //
        f = sample + ".snps.tsv.gz";
        gz_fh = gzopen(f.c_str(), "rb");
        if (!gz_fh) {
            cerr << " Unable to open '" << sample << "'\n";
            return 0;
        }
        #if ZLIB_VERNUM >= 0x1240
        gzbuffer(gz_fh, libz_buffer_size);
        #endif
        gzip = true;
    }
    cerr << "  Parsing " << f.c_str() << "\n";

    line_num = 1;
    while (fh_status) {
        fh_status = (gzip == true) ? read_gzip_line(gz_fh, &line, &size) : read_line(fh, &line, &size);

        if (!fh_status && strlen(line) == 0)
            continue;
        if (is_comment(line)) continue;

        parse_tsv(line, parts);

        if (parts.size() != num_snps_fields && parts.size() != num_snps_fields - 2) {
            cerr << "Error parsing " << f.c_str() << " at line: " << line_num << ". (" << parts.size() << " fields).\n";
            return 0;
        }

        samp_id = atoi(parts[1].c_str());
        id      = atoi(parts[2].c_str());

        snp         = new SNP;
        snp->col    = atoi(parts[3].c_str());

        if (parts[4] == "O")
            snp->type = snp_type_hom;
        else if (parts[4] == "E")
            snp->type = snp_type_het;
        else
            snp->type = snp_type_unk;

        snp->lratio = atof(parts[5].c_str());
        snp->rank_1 = parts[6].at(0);
        snp->rank_2 = parts[7].at(0) == '-' ? 0 : parts[7].at(0);

        if (parts.size() == 10) {
            if (parts[8].length() == 0 || parts[8].at(0) == '-')
                snp->rank_3 = 0;
            else
                snp->rank_3 = parts[8].at(0);
            if (parts[9].length() == 0 || parts[9].at(0) == '-')
                snp->rank_4 = 0;
            else
                snp->rank_4 = parts[9].at(0);
        }

        if (snpres.count(id) == 0) {
            snpr = new SNPRes(samp_id, id);
            snpres[id] = snpr;
        }
        snpres[id]->snps.push_back(snp);

        line_num++;
    }

    if (gzip)
        gzclose(gz_fh);
    else
        fh.close();

    delete [] line;

    return 1;
}
