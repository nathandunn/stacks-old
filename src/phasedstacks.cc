// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2013, Julian Catchen <jcatchen@uoregon.edu>
//
// This file is part of Stacks.
//
// Stacks is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Stacks is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Stacks.  If not, see <http://www.gnu.org/licenses/>.
//

//
// phasedstacks -- analyse phased data, descended from a Stacks analysis.
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//

#include "phasedstacks.h"

// Global variables to hold command-line options.
file_type in_file_type = unknown;
int       num_threads  = 1;
string    in_path;
string    out_path;
string    out_file;
double    minor_freq_lim = 0.2;

int main (int argc, char* argv[]) {

    parse_command_line(argc, argv);

    cerr << "Looking for ";
    switch(in_file_type) {
    case phase:
    default:
	cerr << "Phase";
	break;
    }
    cerr << " input files.\n";

    //
    // Set the number of OpenMP parallel threads to execute.
    //
    #ifdef _OPENMP
    omp_set_num_threads(num_threads);
    #endif

    vector<pair<int, string> > files;
    if (!build_file_list(files))
	exit(1);

    cerr << "Identified " << files.size() << " files.\n";

    //
    // Open the log file.
    //
    stringstream log;
    log << "phasedstacks.log";
    string log_path = in_path + log.str();
    ofstream log_fh(log_path.c_str(), ofstream::out);

    if (log_fh.fail()) {
        cerr << "Error opening log file '" << log_path << "'\n";
	exit(1);
    }

    for (uint i = 0; i < files.size(); i++) {

	if (files[i].second != "batch_1.groupV.phase") continue;

	PhasedSummary *psum;

	if ((psum = parse_phase(in_path + files[i].second)) == NULL) {
	    cerr << "Unable to parse input files.\n";
	    exit(1);
	}

	//
	// Summarize the genotypes in the populations.
	//
	summarize_phased_genotypes(psum);

	// for (uint j = 0; j < psum->size; j++) {
	//     cerr << "BP: " << psum->nucs[j].bp << "\t"
	// 	 << "A: "  << std::setw(3) << psum->nucs[j].nuc[0] << " "
	// 	 << "C: "  << std::setw(3) << psum->nucs[j].nuc[1] << " "
	// 	 << "G: "  << std::setw(3) << psum->nucs[j].nuc[2] << " "
	// 	 << "T: "  << std::setw(3) << psum->nucs[j].nuc[3] << "\n";
	// }

	//
	// Calculate D'
	//
	cerr << "Calculating D'...";
	calc_dprime(psum);
	cerr << "done.\n";

	write_dprime(in_path + files[i].second, psum);

	//
	// Free the Samples objects
	//
	delete psum;
    }

    log_fh.close();

    return 0;
}

int
calc_dprime(PhasedSummary *psum)
{
    #pragma omp parallel
    {
	char  allele_A, allele_a, allele_B, allele_b;
	float freq_A,     freq_a,   freq_B,   freq_b;
	float freq_AB,   freq_Ab,  freq_aB,  freq_ab;
	float D, min;
	float tot = psum->sample_cnt  * 2.0;

	#pragma omp for schedule(dynamic, 1)
	for (uint i = 0; i < psum->size; i++) {
	    //
	    // Assign nucleotides to allele A, and a.
	    //
	    assign_alleles(psum->nucs[i], allele_A, allele_a, freq_A, freq_a);

	    for (uint j = i+1; j < psum->size; j++) {
		//
		// Assign nucleotides to allele B, and b.
		//
		assign_alleles(psum->nucs[j], allele_B, allele_b, freq_B, freq_b);

		freq_AB = 0.0;
		freq_Ab = 0.0;
		freq_aB = 0.0;
		freq_ab = 0.0;
		D       = 0.0;

		//
		// Tally up haplotype frequencies.
		//
		for (uint k = 0; k < psum->sample_cnt; k++) {

		    if (psum->samples[k].nucs_1[i] == allele_A &&
			psum->samples[k].nucs_1[j] == allele_B)
			freq_AB++;
		    else if (psum->samples[k].nucs_1[i] == allele_A &&
			     psum->samples[k].nucs_1[j] == allele_b)
			freq_Ab++;
		    else if (psum->samples[k].nucs_1[i] == allele_a &&
			     psum->samples[k].nucs_1[j] == allele_B)
			freq_aB++;
		    else if (psum->samples[k].nucs_1[i] == allele_a &&
			     psum->samples[k].nucs_1[j] == allele_b)
			freq_ab++;

		    if (psum->samples[k].nucs_2[i] == allele_A &&
			psum->samples[k].nucs_2[j] == allele_B)
			freq_AB++;
		    else if (psum->samples[k].nucs_2[i] == allele_A &&
			     psum->samples[k].nucs_2[j] == allele_b)
			freq_Ab++;
		    else if (psum->samples[k].nucs_2[i] == allele_a &&
			     psum->samples[k].nucs_2[j] == allele_B)
			freq_aB++;
		    else if (psum->samples[k].nucs_2[i] == allele_a &&
			     psum->samples[k].nucs_2[j] == allele_b)
			freq_ab++;
		}

		freq_AB = freq_AB / tot;
		freq_Ab = freq_Ab / tot;
		freq_aB = freq_aB / tot;
		freq_ab = freq_ab / tot;

		D = freq_AB - (freq_A * freq_B);
		cerr << "D_AB: " << D << "; ";
		D = freq_Ab - (freq_A * freq_b);
		cerr << "D_Ab: " << D << "; ";
		D = freq_aB - (freq_a * freq_B);
		cerr << "D_aB: " << D << "; ";
		D = freq_ab - (freq_a * freq_b);
		cerr << "D_ab: " << D << "\n";
		cerr << "    freq_AB: " << freq_AB << "; freq_Ab: " << freq_Ab << "; freq_aB: " << freq_aB << "; freq_ab: " << freq_ab << "\n";

		if (D > 0) {
		    min = (freq_A * freq_b) < (freq_a * freq_B) ? (freq_A * freq_b) : (freq_a * freq_B);
		    D   = D / min;
		} else {
		    min = (freq_A * freq_B) < (freq_a * freq_b) ? (freq_A * freq_B) : (freq_a * freq_b);
		    D   = (-1 * D) / min;
		}

		psum->dprime[i][j] = D;
	    }
	}
    }

    return 0;
}

int
assign_alleles(NucSum nsum, char &p_allele, char &q_allele, float &p_freq, float &q_freq)
{
    p_allele = 0;
    q_allele = 0;

    uint  i   = 0;
    float tot = 0;

    while (p_allele == 0 && i < 4) {
	if (nsum.nuc[i] > 0) {
	    tot += nsum.nuc[i];
	    switch(i) {
	    case 0:
		p_allele = 'A';
		p_freq   = nsum.nuc[0];
		break;
	    case 1:
		p_allele = 'C';
		p_freq   = nsum.nuc[1];
		break;
	    case 2:
		p_allele = 'G';
		p_freq   = nsum.nuc[2];
		break;
	    case 3:
		p_allele = 'T';
		p_freq   = nsum.nuc[3];
		break;
	    }
	}
	i++;
    }
    while (q_allele == 0 && i < 4) {
	if (nsum.nuc[i] > 0) {
	    tot += nsum.nuc[i];
	    switch(i) {
	    case 1:
		q_allele = 'C';
		q_freq   = nsum.nuc[1];
		break;
	    case 2:
		q_allele = 'G';
		q_freq   = nsum.nuc[2];
		break;
	    case 3:
		q_allele = 'T';
		q_freq   = nsum.nuc[3];
		break;
	    }
	}
	i++;
    }

    p_freq = p_freq / tot;
    q_freq = 1 - p_freq;

    return 0;
}

int
write_dprime(string path, PhasedSummary *psum)
{
    //
    // Write the D' data for plotting as a heatmap.
    //
    string file = path + ".dprime.tsv";

    cerr << "Writing D' data to '" << file << "'\n";

    ofstream fh(file.c_str(), ofstream::out);

    if (fh.fail()) {
        cerr << "Error opening D' file '" << file << "'\n";
	exit(1);
    }

    for (uint i = 0; i < psum->size; i++) {
	for (uint j = i+1; j < psum->size; j++) {

	    if (psum->nucs[i].freq < minor_freq_lim || 
		psum->nucs[j].freq < minor_freq_lim)
		continue;

	    fh << psum->nucs[i].bp << "\t" 
	       << psum->nucs[j].bp << "\t"
	       << std::setprecision(3) << psum->dprime[i][j] << "\n";
	}
    }

    fh.close();

    return 0;
}

int
summarize_phased_genotypes(PhasedSummary *psum)
{
    //
    // Construct a two dimensional array out of all the nucleotide arrays in the samples.
    //
    char **gtypes = new char *[psum->sample_cnt];

    for (uint i = 0; i < psum->sample_cnt; i++) {
	gtypes[i] = psum->samples[i].nucs_1;
    }

    //
    // Sum up the occurences of each nucleotide.
    //
    for (uint i = 0; i < psum->size; i++) {
	for (uint j = 0; j < psum->sample_cnt; j++) {
	    switch(gtypes[j][i]) {
	    case 'A':
		psum->nucs[i].nuc[0]++;
		break;
	    case 'C':
		psum->nucs[i].nuc[1]++;
		break;
	    case 'G':
		psum->nucs[i].nuc[2]++;
		break;
	    case 'T':
		psum->nucs[i].nuc[3]++;
		break;
	    case 'N':
	    default:
		break;
	    }
	}
    }

    //
    // Repeat for the second set of phased genotypes.
    //
    for (uint i = 0; i < psum->sample_cnt; i++) {
	gtypes[i] = psum->samples[i].nucs_2;
    }

    //
    // Sum up the occurences of each nucleotide.
    //
    for (uint i = 0; i < psum->size; i++) {
	for (uint j = 0; j < psum->sample_cnt; j++) {
	    switch(gtypes[j][i]) {
	    case 'A':
		psum->nucs[i].nuc[0]++;
		break;
	    case 'C':
		psum->nucs[i].nuc[1]++;
		break;
	    case 'G':
		psum->nucs[i].nuc[2]++;
		break;
	    case 'T':
		psum->nucs[i].nuc[3]++;
		break;
	    case 'N':
	    default:
		break;
	    }
	}

	//
	// Calculate minor allele frequency.
	//
	float tot  = (float) psum->sample_cnt * 2.0;
	float freq = 0.0;
	for (uint j = 0; j < 4; j++) {
	    if (psum->nucs[i].nuc[j] > 0) {
		freq = (float) psum->nucs[i].nuc[j] / tot;
		psum->nucs[i].freq = freq < psum->nucs[i].freq ? freq : psum->nucs[i].freq;
	    }
	}
    }

    delete [] gtypes;

    return 0;
}

//
// Code to parse fastPhase format. 
//
PhasedSummary * 
parse_phase(string path) 
{
    ifstream    fh;
    char        line[max_len];
    string      buf, filepath;
    const char *p, *q, *end;
    int         i, sindex;

    memset(line, '\0', max_len);

    //
    // Read in the original PHASE export from Stacks to obtain the original base pair positions.
    //
    //
    // Open the file for reading
    //
    filepath = path + ".inp";
    fh.open(filepath.c_str(), ifstream::in);

    if (fh.fail()) {
        cerr << "Error opening input file '" << path << "'\n";
	return NULL;
    }

    cerr << "Parsing " << filepath << "...\n";

    int  num_samples, num_genotypes;
    char bp[id_len];

    //
    // Get the number of samples in the dataset.
    //
    fh.getline(line, max_len);
    num_samples = is_integer(line);

    if (num_samples < 0) {
	cerr << "Unable to find the number of samples, should be the first line.\n";
	return NULL;
    }

    //
    // Get the number of genotypes in the dataset.
    //
    fh.getline(line, max_len);
    num_genotypes = is_integer(line);

    if (num_genotypes < 0) {
	cerr << "Unable to find the number of genotypes, should be the second line.\n";
	return NULL;
    }

    PhasedSummary *psum = new PhasedSummary(num_samples, num_genotypes);

    //
    // Get the set of base pair positions.
    //
    buf.clear();
    do {
	fh.clear();
	fh.getline(line, max_len);
	buf += line;
    } while (fh.fail() && !fh.bad() && !fh.eof());

    i   = 0;
    p   = buf.c_str();
    end = p + buf.length();

    if (*p != 'P') {
	cerr << "Unable to locate line of basepair positions, should be the third line.\n";
	delete psum;
	return NULL;
    }
    for (p += 2, q = p; p < end; p++, q++) {
	while (*q != ' ' && q < end) {
	    q++; 
	}
	strncpy(bp, p, q - p);
	bp[q - p] = '\0';
	psum->nucs[i].bp = is_integer(bp);

	if (psum->nucs[i].bp < 0) {
	    cerr << "Unable to parse base pair positions.\n";
	    delete psum;
	    return NULL;
	}

	i++;
	p = q;
    }

    fh.close();

    //
    // Open the file for reading
    //
    filepath = path + "_hapguess_switch.out";
    fh.open(filepath.c_str(), ifstream::in);

    if (fh.fail()) {
        cerr << "Error opening input file '" << path << "'\n";
	return NULL;
    }

    cerr << "Parsing " << filepath << "...\n";

    //
    // Read from the "*_hapguess_switch.out" file until we hit the genotypes section
    // marked by the string "BEGIN GENOTYPES".
    //
    do {
        fh.getline(line, max_len);

        if (!fh.good()) {
	    cerr << "Unable to find file section entitled 'BEGIN GENOTYPES'\n";
	    delete psum;
            return NULL;
	}

    } while (strcmp(line, "BEGIN GENOTYPES") != 0);

    //
    // Now read lines from the file in groups of three:
    //   1. Sample label
    //   2. Phased genotypes from chromosome 1
    //   3. Phased genotypes from chromosome 2
    // Stop reading individuals when we encounter the string, "END GENOTYPES".
    //
    fh.getline(line, max_len);

    do {
	//
	// Create a new Sample object and store the sample label.
	//
	sindex = psum->add_sample(line);

	//
	// Get the first set of phased genotypes.
	//
	buf.clear();
	do {
	    fh.clear();
	    fh.getline(line, max_len);
	    buf += line;
	} while (fh.fail() && !fh.bad() && !fh.eof());

	//
	// Count the number of genotypes on this line (they should be space deliniated).
	//
	i = 0;
	for (p = buf.c_str(); *p != '\0'; p++)
	    if (*p != ' ') psum->samples[sindex].size++;
	//
	// Store the genotypes into our internal buffer.
	//
	psum->samples[sindex].nucs_1 = new char[psum->samples[sindex].size];
	for (p = buf.c_str(); *p != '\0'; p++) {
	    if (*p == ' ') continue;
	    psum->samples[sindex].nucs_1[i] = *p;
	    i++;
	}

	// len = strlen(line);
	// if (line[len - 1] == '\r') line[len - 1] = '\0';

	//
	// Get the second set of phased genotypes.
	//
	buf.clear();
	do {
	    fh.clear();
	    fh.getline(line, max_len);
	    buf += line;
	} while (fh.fail() && !fh.bad() && !fh.eof());

	i = 0;
	psum->samples[sindex].nucs_2 = new char[psum->samples[sindex].size];
	for (p = buf.c_str(); *p != '\0'; p++) {
	    if (*p == ' ') continue;
	    psum->samples[sindex].nucs_2[i] = *p;
	    i++;
	}

	//
	// Get the sample label of the next record.
	//
	fh.getline(line, max_len);

    } while (strcmp(line, "END GENOTYPES") != 0 && fh.good());

    fh.close();

    return psum;
}

int build_file_list(vector<pair<int, string> > &files) {
    vector<string> parts;
    string pattern;

    //
    // Read all the files from the Stacks directory.
    //
    uint   pos;
    string file;
    struct dirent *direntry;

    DIR *dir = opendir(in_path.c_str());

    if (dir == NULL) {
	cerr << "Unable to open directory '" << in_path << "' for reading.\n";
	exit(1);
    }

    switch(in_file_type) {
    case phase:
    default:
	pattern = "_hapguess_switch.out";
	break;
    }

    while ((direntry = readdir(dir)) != NULL) {
	file = direntry->d_name;

	if (file == "." || file == "..")
	    continue;

	pos = file.rfind(pattern);
	if (pos < file.length())
	    files.push_back(make_pair(1, file.substr(0, pos)));
    }

    closedir(dir);

    if (files.size() == 0) {
	cerr << "Unable to locate any input files to process within '" << in_path << "'\n";
	return 0;
    }

    return 1;
}

int parse_command_line(int argc, char* argv[]) {
    int c;
     
    while (1) {
	static struct option long_options[] = {
	    {"help",        no_argument,       NULL, 'h'},
            {"version",     no_argument,       NULL, 'v'},
            {"infile_type", required_argument, NULL, 't'},
	    {"num_threads", required_argument, NULL, 'p'},
	    {"in_path",     required_argument, NULL, 'P'},
	    {0, 0, 0, 0}
	};
	
	// getopt_long stores the option index here.
	int option_index = 0;
     
	c = getopt_long(argc, argv, "hvAt:P:p:", long_options, &option_index);
     
	// Detect the end of the options.
	if (c == -1)
	    break;
     
	switch (c) {
	case 'h':
	    help();
	    break;
	case 'p':
	    num_threads = atoi(optarg);
	    break;
	case 'P':
	    in_path = optarg;
	    break;
     	case 't':
            if (strcasecmp(optarg, "phase") == 0)
                in_file_type = phase;
            else
                in_file_type = unknown;
	    break;
        case 'v':
            version();
            break;
	case '?':
	    // getopt_long already printed an error message.
	    help();
	    break;
	default:
	    help();
	    abort();
	}
    }

    if (in_path.length() == 0) {
	cerr << "You must specify a path to the directory containing Stacks output files.\n";
	help();
    }

    if (in_path.at(in_path.length() - 1) != '/') 
	in_path += "/";

    return 0;
}

void version() {
    std::cerr << "phasedstacks " << VERSION << "\n\n";

    exit(0);
}

void help() {
    std::cerr << "phasedstacks " << VERSION << "\n"
              << "phasedstacks -P path -t file_type [-p threads] [-v] [-h]" << "\n"
	      << "  P: path to the phased Stacks output files.\n"
	      << "  t: input file type. Supported types: phase.\n"
	      << "  p: number of processes to run in parallel sections of code.\n"
	      << "  v: print program version." << "\n"
	      << "  h: display this help messsage." << "\n\n";

    exit(0);
}
