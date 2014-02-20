// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright (c) 2014 University of Oregon
// Created by Julian Catchen <jcatchen@uoregon.edu>
//

#ifndef __GZFASTQ_H__
#define __GZFASTQ_H__

#ifdef HAVE_LIBZ

#include <errno.h>
#include <zlib.h>
#include "input.h"

class GzFastq: public Input {
    gzFile gz_fh;

public:
    GzFastq(string path) : Input() { 
	this->gz_fh = gzopen(path.c_str(), "rb");
	if (!this->gz_fh) {
	    cerr << "Failed to open gzipped file '" << path << "': " << strerror(errno) << ".\n";
	    exit(EXIT_FAILURE);
	}
    };
    GzFastq(const char *path) : Input() { 
	this->gz_fh = gzopen(path, "rb");
	if (!this->gz_fh) {
	    cerr << "Failed to open gzipped file '" << path << "': " << strerror(errno) << ".\n";
	    exit(EXIT_FAILURE);
	}
    };
    ~GzFastq() {
	gzclose(this->gz_fh);
    };
    Seq *next_seq();
    int  next_seq(Seq &s);
};

Seq *GzFastq::next_seq() {
    char *res = NULL;

    //
    // Check the contents of the line buffer. When we finish reading a FASTQ record
    // the buffer will either contain whitespace or the header of the next FASTQ
    // record.
    //
    this->line[0] = '\0';
    do {
	res = gzgets(this->gz_fh, this->line, max_len);
    } while (this->line[0] != '@' && res != NULL);

    if (res == NULL) {
	return NULL;
    }

    //
    // Check if there is a carraige return in the buffer
    //
    uint len = strlen(this->line);
    if (len > 0 && this->line[len - 1] == '\n') this->line[len - 1] = '\0';
    if (len > 0 && this->line[len - 2] == '\r') this->line[len - 2] = '\0';

    //
    // Initialize the Seq structure and store the FASTQ ID
    //
    Seq *s = new Seq;
    s->id = new char[strlen(this->line) + 1];
    strcpy(s->id, this->line + 1);

    //
    // Read the sequence from the file
    //
    gzgets(this->gz_fh, this->line, max_len);

    if (gzeof(this->gz_fh)) {
	return NULL;
    }

    len = strlen(this->line);
    if (len > 0 && this->line[len - 1] == '\n') this->line[len - 1] = '\0';
    if (len > 0 && this->line[len - 2] == '\r') this->line[len - 2] = '\0';

    s->seq = new char[len + 1];
    strcpy(s->seq, this->line);

    //
    // Read the repeat of the ID
    //
    this->line[0] = '\0';
    res = gzgets(this->gz_fh, this->line, max_len);

    if (this->line[0] != '+' || res == NULL) {
	return NULL;
    }

    //
    // Read the quality score from the file
    //
    this->line[0] = '\0';
    res = gzgets(this->gz_fh, this->line, max_len);

    if (res == NULL && strlen(this->line) == 0) {
	return NULL;
    }

    len = strlen(this->line);
    if (len > 0 && this->line[len - 1] == '\n') this->line[len - 1] = '\0';
    if (len > 0 && this->line[len - 2] == '\r') this->line[len - 2] = '\0';

    s->qual = new char[len + 1];
    strcpy(s->qual, this->line);

    //
    // Clear the line buffer so it is set up for the next record. If a '@'
    // appears in the quality scores read, it will break parsing next time 
    // it is called.
    //
    this->line[0] = '\0';

    return s;
}

int GzFastq::next_seq(Seq &s) {
    char *res = NULL;

    //
    // Check the contents of the line buffer. When we finish reading a FASTQ record
    // the buffer will either contain whitespace or the header of the next FASTQ
    // record.
    //
    this->line[0] = '\0';
    do {
	res = gzgets(this->gz_fh, this->line, max_len);
    } while (this->line[0] != '@' && res != NULL);

    if (res == NULL) {
	return 0;
    }

    //
    // Check if there is a carraige return in the buffer
    //
    uint len = strlen(this->line);
    if (len > 0 && this->line[len - 1] == '\n') this->line[len - 1] = '\0';
    if (len > 0 && this->line[len - 2] == '\r') this->line[len - 2] = '\0';

    //
    // Store the FASTQ ID
    //
    strcpy(s.id, this->line + 1);

    //
    // Read the sequence from the file
    //
    this->line[0] = '\0';
    res = gzgets(this->gz_fh, this->line, max_len);

    if (res == NULL) {
	return 0;
    }

    len = strlen(this->line);
    if (len > 0 && this->line[len - 1] == '\n') this->line[len - 1] = '\0';
    if (len > 0 && this->line[len - 2] == '\r') this->line[len - 2] = '\0';

    strcpy(s.seq, this->line);

    //
    // Read the repeat of the ID
    //
    this->line[0] = '\0';
    res = gzgets(this->gz_fh, this->line, max_len);

    if (this->line[0] != '+' || res == NULL) {
	return 0;
    }

    //
    // Read the quality score from the file
    //
    this->line[0] = '\0';
    res = gzgets(this->gz_fh, this->line, max_len);

    if (res == NULL && strlen(this->line) == 0) {
	return 0;
    }

    len = strlen(this->line);
    if (len > 0 && this->line[len - 1] == '\n') this->line[len - 1] = '\0';
    if (len > 0 && this->line[len - 2] == '\r') this->line[len - 2] = '\0';

    strcpy(s.qual, this->line);

    //
    // Clear the line buffer so it is set up for the next record. If a '@'
    // appears in the quality scores read, it will break parsing next time 
    // it is called.
    //
    this->line[0] = '\0';

    return 1;
}

#else  // If HAVE_LIBZ is undefined and zlib library is not present.

#include "input.h"

class GzFastq: public Input {
 public:
    GzFastq(const char *path) : Input() { cerr << "Gzip support was not enabled when Stacks was compiled.\n"; };
    GzFastq(string path) : Input() { cerr << "Gzip support was not enabled when Stacks was compiled.\n"; };
    ~GzFastq() {};
    Seq *next_seq()      { return NULL; };
    int  next_seq(Seq &) { return 0; };
};

#endif // HAVE_LIBZ

#endif // __GZFASTQ_H__