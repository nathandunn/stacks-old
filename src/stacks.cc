// -*-mode:c++; c-style:k&r; c-basic-offset:4;-*-
//
// Copyright 2010, Julian Catchen <jcatchen@uoregon.edu>
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
// stacks.cc -- routines for the stack-holding containers
//
// Julian Catchen
// jcatchen@uoregon.edu
// University of Oregon
//
// $Id$
//
#include "stacks.h"

Rem::Rem() { 
    this->id         = 0;
    this->seq        = NULL; 
    this->utilized   = false;
}

Rem::Rem(int id, uint seq_id, DNASeq *seq) { 
    this->id       = id;
    this->utilized = false;

    this->map.push_back(seq_id);

    this->seq = new DNASeq(seq->size, seq->s);
}

int Rem::add_id(uint id) {
    this->map.push_back(id);

    return 0;
}

int Rem::add_seq(const DNASeq *seq) {
    if (this->seq != NULL)
	delete this->seq;

    this->seq = new DNASeq(seq->size, seq->s);

    return 0;
}

int Rem::add_seq(const char *seq) {
    if (this->seq != NULL)
	delete this->seq;

    this->seq = new DNASeq(strlen(seq), seq);

    return 0;
}

int PStack::add_id(const char *id) {
    char *f = new char[strlen(id) + 1];
    strcpy(f, id);
    this->map.push_back(f);

    return 0;
}

int PStack::add_seq(const char *seq) {
    if (this->seq != NULL)
	delete this->seq;

    this->len = strlen(seq);
    this->seq = new DNANSeq(this->len, seq);

    return 0;
}

int PStack::add_seq(DNANSeq *seq) {
    if (this->seq != NULL)
	delete this->seq;

    this->seq = new DNANSeq(seq->size(), seq->s);

    return 0;
}

int Stack::add_id(uint id) {
    this->map.push_back(id);

    return 0;
}

int Stack::add_seq(const char *seq) {
    if (this->seq != NULL)
	delete this->seq;

    this->seq = new DNASeq(strlen(seq), seq);

    return 0;
}

int Stack::add_seq(const DNASeq *seq) {
    if (this->seq != NULL)
	delete this->seq;

    this->seq = new DNASeq(seq->size, seq->s);

    return 0;
}
