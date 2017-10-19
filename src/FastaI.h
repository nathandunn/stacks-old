// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010-2013, Julian Catchen <jcatchen@uoregon.edu>
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

#ifndef __FASTAI_H__
#define __FASTAI_H__

#include "input.h"

class Fasta: public Input {
    string buf;

 public:
    Fasta(const char *path) : Input(path) {
        if (fh.peek() != '>') {
            cerr << "Error: '" << path << "': not in fasta format (expected '>').\n";
            throw exception();
        }
    };
    Fasta(string path) : Fasta(path.c_str()) { };
    ~Fasta() {};
    Seq *next_seq();
    int  next_seq(Seq &);
};

Seq *Fasta::next_seq() {
    //
    // Check the contents of the line buffer. When we finish reading a FASTA record
    // the buffer will either contain whitespace or the header of the next FAST
    // record.
    //
    while (this->line[0] != '>' && this->fh.good() ) {
        this->fh.getline(this->line, max_len);
    }

    if (!this->fh.good()) {
        return NULL;
    }

    //
    // Check if there is a carraige return in the buffer
    //
    uint len = strlen(this->line);
    if (this->line[len - 1] == '\r') this->line[len - 1] = '\0';

    //
    // Initialize the Seq structure and store the FASTA ID
    //
    Seq *s = new Seq;

    //
    // Check if the ID line of the FASTA file has a comment after the ID.
    //
    const char *p, *q;
    p = this->line + 1;
    for (q = this->line; *q != '\0' && *q != ' ' && *q != '\t'; q++);

    if (*q == '\0') {
        // Comment not present.
        s->id = new char[len + 1];
        strcpy(s->id, p);

    } else {
        // Comment present.
        int l = q - p;
        assert(l > 0);
        s->id = new char[l + 1];
        strncpy(s->id, p, l);
        s->id[l] = '\0';

        q++;
        p = q;
        for (; *q != '\0'; q++);
        l = q - p;
        assert(l > 0);
        s->comment = new char[l + 1];
        strncpy(s->comment, p, l);
        s->comment[l] = '\0';
    }

    //
    // Read the sequence from the file -- keep reading lines until we reach the next
    // record or the end of file.
    //
    this->fh.getline(this->line, max_len);

    while (this->line[0] != '>' && this->fh.good()) {
        len = strlen(this->line);
        if (this->line[len - 1] == '\r') this->line[len - 1] = '\0';

        this->buf    += this->line;
        this->line[0] = '\0';
        this->fh.getline(this->line, max_len);
    }

    if (this->fh.eof()) {
        len = strlen(this->line);
        if (this->line[len - 1] == '\r') this->line[len - 1] = '\0';

        this->buf += this->line;
    }

    s->seq = new char[this->buf.length() + 1];
    strcpy(s->seq, this->buf.c_str());
    this->buf.clear();

    return s;
}

int Fasta::next_seq(Seq &s) {
    //
    // Check the contents of the line buffer. When we finish reading a FASTA record
    // the buffer will either contain whitespace or the header of the next FAST
    // record.
    //
    while (this->line[0] != '>' && this->fh.good() ) {
        this->fh.getline(this->line, max_len);
    }

    if (!this->fh.good()) {
        return 0;
    }

    //
    // Check if there is a carraige return in the buffer
    //
    uint len = strlen(this->line);
    if (this->line[len - 1] == '\r') this->line[len - 1] = '\0';

    //
    // Check if the ID line of the FASTA file has a comment after the ID.
    //
    const char *p, *q;
    p = this->line + 1;
    for (q = this->line; *q != '\0' && *q != ' ' && *q != '\t'; q++);

    if (*q == '\0') {
        // Comment not present.
        strncpy(s.id, p, id_len);
        s.id[id_len - 1] = '\0';

    } else {
        // Comment present.
        int l = q - p;
        assert(l > 0);
        strncpy(s.id, p, l);
        s.id[l] = '\0';

        q++;
        p = q;
        for (; *q != '\0'; q++);
        l = q - p;
        assert(l > 0);
        strncpy(s.comment, p, l);
        s.comment[l] = '\0';
    }

    //
    // Read the sequence from the file -- keep reading lines until we reach the next
    // record or the end of file.
    //
    this->fh.getline(this->line, max_len);

    while (this->line[0] != '>' && this->fh.good()) {
        len = strlen(this->line);
        if (len > 0 && this->line[len - 1] == '\r') this->line[len - 1] = '\0';

        this->buf    += this->line;
        this->line[0] = '\0';
        this->fh.getline(this->line, max_len);
    }

    if (this->fh.eof()) {
        len = strlen(this->line);
        if (len > 0 && this->line[len - 1] == '\r') this->line[len - 1] = '\0';

        this->buf += this->line;
    }

    strcpy(s.seq, this->buf.c_str());
    this->buf.clear();

    return 1;
}

#endif // __FASTAI_H__
