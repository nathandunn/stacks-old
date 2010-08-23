// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-

#ifndef __FASTA_H__
#define __FASTA_H__

#include "input.h"

class Fasta: public Input {
    string buf;

 public:
    Fasta(const char *path) : Input(path) {};
    ~Fasta() {};
    Seq *next_seq();
};

Seq *Fasta::next_seq() {
    //
    // Check the contents of the line buffer. When we finish reading a FASTA record
    // the buffer will either contain whitespace or the header of the next FAST
    // record.
    //
    while (line[0] != '>' && this->fh.good() ) {
	this->fh.getline(this->line, max_len);
    }

    if (!this->fh.good()) {
	return NULL;
    }

    //
    // Initialize the Seq structure and store the FASTA ID
    //
    Seq *s = new Seq;
    s->id = new char[strlen(this->line) + 1];
    strcpy(s->id, this->line + 1);

    //
    // We will also load the FASTA ID into the Seq object as a fake genomic location
    // so that the stacks are built correctly in the pstacks program.
    //
    s->chr = new char[strlen(this->line) + 1];
    strcpy(s->chr, this->line + 1);
    s->bp = 0;

    //
    // Read the sequence from the file -- keep reading lines until we reach the next
    // record or the end of file.
    //
    this->fh.getline(this->line, max_len);

    while (this->line[0] != '>' && this->fh.good()) {
	this->buf += this->line;
	this->fh.getline(this->line, max_len);
    }

    if (this->fh.eof()) {
	this->buf += this->line;
    }

    s->seq = new char[this->buf.length() + 1];
    strcpy(s->seq, this->buf.c_str());
    this->buf.clear();

    return s;
}

#endif // __FASTA_H__
