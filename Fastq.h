// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-

#ifndef __FASTQ_H__
#define __FASTQ_H__

#include "input.h"

class Fastq: public Input {

public:
    Fastq(const char *path) : Input(path) {};
    ~Fastq() {};
    Seq *next_seq();
};

Seq *Fastq::next_seq() {
    //
    // Check the contents of the line buffer. When we finish reading a FASTQ record
    // the buffer will either contain whitespace or the header of the next FASTQ
    // record.
    //
    while (line[0] != '@' && this->fh.good() ) {
	this->fh.getline(this->line, max_len);
    }

    if (!this->fh.good()) {
	return NULL;
    }

    //
    // Initialize the Seq structure and store the FASTQ ID
    //
    Seq *s = new Seq;
    s->id = new char[strlen(this->line) + 1];
    strcpy(s->id, this->line + 1);

    //
    // We will also load the FASTQ ID into the Seq object as a fake genomic location
    // so that the stacks are built correctly in the pstacks program.
    //
    s->chr = new char[strlen(this->line) + 1];
    strcpy(s->chr, this->line + 1);
    s->bp = 0;

    //
    // Read the sequence from the file
    //
    this->fh.getline(this->line, max_len);

    if (!this->fh.good()) {
	return NULL;
    }

    s->seq = new char[strlen(this->line) + 1];
    strcpy(s->seq, this->line);

    //
    // Read the repeat of the ID
    //
    this->fh.getline(this->line, max_len);

    if (this->line[0] != '+' || !this->fh.good()) {
	return NULL;
    }

    //
    // Read the quality score from the file
    //
    this->fh.getline(this->line, max_len);

    if (!this->fh.good() && !this->fh.eof()) {
	return NULL;
    }

    s->qual = new char[strlen(this->line) + 1];
    strcpy(s->qual, this->line);

    return s;
}

#endif // __FASTQ_H__
