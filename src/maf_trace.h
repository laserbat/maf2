/*
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

/* Proof tracing for automata and kbprog.

   When switched on (the "-trace file" option) every equation MAF stores,
   in the tree or in the pool, is written to the trace file together with a
   derivation: a list of rewriting steps that turns its left hand side into
   its right hand side, each step citing an equation written out earlier.
   The rules of the last rewriting system MAF writes out are listed at the
   end, so that a proof of any of them can be reconstructed by following the
   citations back to the axioms.

   The format is the one KBMAG's kbprog writes under its own -trace option
   (see lib/kbtrace.h in https://github.com/laserbat/kbmag, branch
   proof-trace), so the same tools read both:

     #KBTRACE 1                   the version of this format
     #ordering <n>                the word ordering, numbered as in KBMAG
                                  (0 shortlex, 1 recursive, 2 rt_recursive,
                                   3 wtlex, 4 wreathprod, 5 anything else)
     #resumed <0|1>               whether the run resumed an earlier one
     #gens <n> <name> ...         the generators, in order
     #inverses <n> <g> ...        the inverse of each generator, 0 for none
     E <pid> <lhs> <rhs> <origin> : <step> ...
     #confluent <0|1>
     #lostinfo <0|1>
     F <n> <pid> ...              the rules of the rewriting system written out

   A word is a comma separated list of generator numbers, counting from 1,
   or "." for the empty word.  A step is "<pid>@<pos><dir>": rewrite the
   subword at the 0-based position <pos> using equation <pid>, from its left
   hand side to its right hand side if <dir> is ">" and the other way round
   if it is "<".  Applying an E line's steps to its <lhs> in turn gives its
   <rhs>, and a step always cites an equation written out earlier.  A step
   citing equation 0 is a hole: a part of the derivation that MAF did not
   record in a form that can be turned into rewriting steps.

   Equations are identified by a persistent id (pid).  A node of MAF's tree
   that is the left hand side of an equation names different equations over
   the run, as its right hand side is replaced, so the trace keeps a map from
   nodes to the pid of their current equation, and resolves the equation a
   rewrite used at the moment the rewrite is recorded.  Since a pid names a
   true equation for good, any pid with the right words will do as a
   citation: equations in the pool, and those the node map has lost track
   of, are found by their words.

   Each derivation is checked as it is written: the steps are replayed on
   the words, and a derivation that does not replay is written as a hole
   instead.  The trace therefore never claims a step that does not hold.
*/

#pragma once
#ifndef MAF_TRACE_INCLUDED
#define MAF_TRACE_INCLUDED 1

#include <stdio.h>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#ifndef MAFBASE_INCLUDED
#include "mafbase.h"
#endif
#ifndef NODEBASE_INCLUDED
#include "nodebase.h"
#endif

class MAF;
class Node_Manager;
class Rewriter_Machine;
class Rewriting_System;
class Working_Equation;
class Full_Derivation;
class Word;

struct Trace_Step
{
  int pid;
  int pos;
  int dir; /* 1 = lhs to rhs, -1 = rhs to lhs */
};

typedef std::vector<Trace_Step> Trace_Chain;
typedef std::vector<Ordinal> Trace_Word;

class Proof_Trace
{
  private:
    FILE * file;
    std::string filename;
    int next_pid;
    unsigned long gaps;
    bool debug;
    Ordinal nr_generators;
    std::vector<Ordinal> inverse_of;
    /* pid -> the equation's words, for replaying derivations */
    std::vector<Trace_Word> lhs_of;
    std::vector<Trace_Word> rhs_of;
    /* lhs and rhs -> the latest pid with those words */
    std::unordered_map<std::string,int> by_content;
    /* the equations that may be written as axioms */
    std::unordered_set<std::string> given;
    std::unordered_map<Node_ID,int> node_pid;
    /* right_inverse[g] is the pid of some equation g*w = 1,
       left_inverse[g] of some equation w*g = 1, or 0 */
    std::vector<int> right_inverse;
    std::vector<int> left_inverse;
    /* the rules of the last rewriting system written out */
    std::vector<int> final_pids;
    bool have_final;
    bool final_confluent;
    bool final_lostinfo;
    bool resumed;
    bool closed;
  public:
    Proof_Trace(const char * filename);
    ~Proof_Trace();
    static void set_filename(const char * filename);
    static const char * wanted();
    /* Opens the file and writes the header. Returns false, having said why,
       if tracing cannot be done for this MAF object */
    bool start(MAF & maf);
    bool started() const
    {
      return file != 0 || closed;
    }
    void close();
    /* Called by Presentation::add_axiom() for each axiom it is given, so
       that the trace can tell these from equations MAF itself re-adds as
       axioms */
    void given_axiom(const Word & lhs,const Word & rhs);
    void set_resumed()
    {
      resumed = true;
    }

    /* The pid of the equation that is currently at the equation node n,
       with its raw right hand side, or 0 if there is none we know of. */
    int node_equation(const Node_Manager & nm,Node_ID n);
    /* Append the steps that rewrite lhs(n) at pos to the fully reduced
       right hand side of n (following right hand sides that are themselves
       left hand sides). Steps that cannot be resolved get pid 0. */
    void node_chain(const Node_Manager & nm,Node_ID n,int pos,Trace_Chain * chain);
    /* The latest pid of an equation with these words, or 0 */
    int equation(const Word & lhs,const Word & rhs) const;

    int axiom(const Word & lhs,const Word & rhs);
    /* Write an equation with the given derivation, after checking that
       the derivation replays. Returns its pid. */
    int emit(const char * origin,const Word & lhs,const Word & rhs,
             const Trace_Chain & chain);
    /* The derivation of a Working_Equation: set up its starting point,
       and later turn everything it recorded into a chain */
    void begin(Working_Equation & we);
    int emit_equation(Working_Equation & we);

    void set_node(Node_ID n,int pid)
    {
      node_pid[n] = pid;
    }
    void final_system(const Rewriting_System & rws,const Node_Manager & nm,
                      bool lostinfo);
  private:
    void put_word(const Ordinal * values,size_t length);
    void put_word(const Trace_Word & w)
    {
      put_word(w.empty() ? 0 : &w[0],w.size());
    }
    int write(const char * origin,const Trace_Word & lhs,const Trace_Word & rhs,
              const Trace_Chain & chain);
    static std::string key(const Trace_Word & lhs,const Trace_Word & rhs);
    int content(const Trace_Word & lhs,const Trace_Word & rhs) const;
    bool replay(Trace_Word w,const Trace_Chain & chain,const Trace_Word & target,
                std::string * why) const;
    bool apply(Trace_Word * w,const Trace_Step & step,std::string * why) const;
    bool cancel(Trace_Word * lhs,Trace_Word * rhs,Trace_Chain * chain,
                size_t leading,size_t trailing,std::string * why) const;
    void complain(const char * what,const std::string & why) const;
    bool justify(const Trace_Word & lhs,const Trace_Word & rhs,
                 const std::vector<Trace_Chain> & candidates,
                 Trace_Chain * chain) const;
    void conjugates(const Node_Manager & nm,const Full_Derivation & d,
                    const Trace_Word & l,const Trace_Word & r,
                    std::vector<Trace_Chain> * candidates,int * p1,int * p2);
};

/* Non-zero only while a trace is being written */
extern Proof_Trace * proof_trace;

#endif
