/*
  Copyright 2008,2009,2010 Alun Williams
  This file is part of MAF.
  MAF is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  MAF is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with MAF.  If not, see <http://www.gnu.org/licenses/>.
*/


/*
$Log: maf_rws.h $
Revision 1.8  2010/06/10 13:58:10Z  Alun
All tabs removed again
Revision 1.7  2010/05/16 23:24:19Z  Alun
RWSC_NO_H flag added. is_confluent() method added. Long comment changed
Revision 1.6  2009/09/13 10:21:08Z  Alun
printf() style functions cleaned up considerably
code now compiles cleanly in g++ both 32-bit and 64-bit
all tabs removed from source files
Revision 1.5  2008/11/04 22:27:50Z  Alun
Re-check in of source from web site after catastrophic loss of correct state
of RCS archive.
Revision 1.6  2008/11/04 23:27:49Z  Alun
Added facility to use RWS constructor as way of building word-acceptor
Revision 1.5  2008/10/30 09:25:24Z  Alun
Reworked to use new style packed words
Revision 1.4  2008/09/29 16:09:21Z  Alun
Switch to using Hash+Collection_Manager as replacement for Indexer.
Currently this is about 10% slower, but this is a more flexible
architecture.
Revision 1.3  2007/11/15 22:58:12Z  Alun
*/
#pragma once
#ifndef MAF_RWS_INCLUDED
#define MAF_RWS_INCLUDED 1

//Classes defined in this header file
class Rewriting_System;

//Classes referred to and defined elsewhere
class Rewriter_Machine;
class Hash;
class MAF;

#ifndef FSA_INCLUDED
#include "fsa.h"
#endif

#ifndef ARRAYBOX_INCLUDED
#include "arraybox.h"
#endif

#ifndef MAFWORD_INCLUDED
#include "mafword.h"
#endif

/* This class is forced into a slightly unsatisfactory architecture in order
   to support KBMAG format rewriting system+index automata. KBMAG stores the
   rewriting system and the index automaton in separate files: a .kbprog file
   containing equations and a .reduce file containing the FSA. This makes
   reading and writing Rewriting_System messy compared with an ordinary FSA.
   Another annoying feature is that KBMAG wants the whole of the LHS
   present for each equation, even though for reduction all we need to know
   is how long it is. A format that stored just the RHS with padding
   characters to indicate the LHS length would be elegant and could have been
   stored using a labelled FSA structure. The LHS could always be recovered
   by seeing what we have read. On the other hand the KBMAG format renders
   the .reduce file essentially superfluous since it can be recreated from
   the .kbprog file alone, and the .kbprog file would be needed if one wanted
   to do reduction some other way, for example using the Rabin-Karp algorithm.
   But, for some presentations the rewriting system would have fewer
   states, and fewer equations if we could omit the LHS.

   We also use this class to generate .kbprog files alone from a MAF object.
   In that case the create_fsa parameter passed to the constructor is false.
*/

const unsigned RWSC_MINIMAL = 1;
const unsigned RWSC_NEED_FSA = 2;
const unsigned RWSC_CRITICAL_ONLY = 4;
const unsigned RWSC_FSA_ONLY = 8;
const unsigned RWSC_NO_H = 16;

class Rewriting_System : public Delegated_FSA
{
  private:
    Hash & eqn;
    Array_Of_Data packed_lhses;
    Array_Of_Data packed_rhses;
    bool confluent;
    Ordinal *inverses;
    /* the machine a minimal system (what goes in the .kbprog file) was
       built from, or 0: a proof trace needs it to resolve the equations
       print() writes out */
    Rewriter_Machine * source;
  public:
    // if the taker parameter is non-zero, which only makes sense when
    // the RWSC_NEED_FSA flag is set, the FSA created by the constructor
    // is given to taker, and it becomes the caller's responsibility to
    // deal with it. In this case the Rewriting_System is left without an
    // FSA, so it cannot perform reduction.
    // But this is useful because it saves having to duplicate the code
    // for construction, and it also saves us having to duplicate the
    // constructed FSA at run time
    Rewriting_System(Rewriter_Machine * rm,
                     unsigned flags = RWSC_MINIMAL|RWSC_NEED_FSA,
                     FSA_Simple **taker = 0);
    Rewriting_System(const MAF & maf,String fsa_filename = 0);
    Rewriting_System(const Rewriting_System & other);
    virtual ~Rewriting_System();
    void print(Output_Stream * kb_stream,Output_Stream *rws_stream) const;
    void set_name(String name_);
    static Rewriting_System * create(String kb_filename,String fsa_filename,
                                     Container * container,
                                     bool must_succeed = true);
#if defined(HASH_INCLUDED) && defined(MAFWORD_INCLUDED)
    State_Count equation_count() const
    {
      /*The return value is actually 1 more than the number of equations:
        If there are n equations they are numbered from 1..n, and the return
        value is n+1 */
      return eqn.count();
    }
    bool is_confluent() const
    {
      return confluent;
    }
    const Byte * packed_lhs(Element_ID eqn_nr) const
    {
      return (const Byte *) packed_lhses.get(eqn_nr);
    }
    const Byte * packed_rhs(Element_ID eqn_nr) const
    {
      return (const Byte *) packed_rhses.get(eqn_nr);
    }
    Word_Length lhs_length(Element_ID eqn_nr) const
    {
      return base_alphabet.extra_packed_word_length(packed_lhs(eqn_nr));
    }
    Word_Length read_lhs(Word * lhs_word,Element_ID eqn_nr) const
    {
      return lhs_word->extra_unpack(packed_lhs(eqn_nr));
    }
    Word_Length read_rhs(Word * rhs_word,Element_ID eqn_nr) const
    {
      return rhs_word->extra_unpack(packed_rhs(eqn_nr));
    }
#endif
  private:
    /* declared to show that print() above is intentional */
    void print(Output_Stream *,unsigned) const {};
};

#ifdef MAF_INCLUDED
class RWS_Reducer : public Word_Reducer
{
  private:
    State_ID *state;
    Word_Length max_stack;
    Ordinal_Word stack_word;
  public:
    Rewriting_System & rws;
    RWS_Reducer(Rewriting_System &rws_) :
      rws(rws_),
      state(0),
      stack_word(rws_.base_alphabet)
    {
      max_stack = 0;
    }
    ~RWS_Reducer()
    {
      if (state)
        delete [] state;
    }
    unsigned reduce(Word * word,const Word & start_word,
                    unsigned flags = 0,const FSA * wa = 0);
  private:
    unsigned hard_reduce(Word * word,unsigned flags,Word_Length valid_length,Ordinal_Word & rhs_word);
};
#endif

#endif
