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
$Log: maf_we.h $
Revision 1.16  2010/04/10 19:05:30Z  Alun
Lots of changes needed for new style Node_Manager interface
Revision 1.15  2009/09/16 07:06:00Z  Alun
Additional source code changes needed to get clean compilation on latest GNU

Revision 1.14  2009/09/12 20:03:39Z  Alun
printf() style functions cleaned up considerably
code now compiles cleanly in g++ both 32-bit and 64-bit
all tabs removed from source files
Revision 1.13  2008/11/26 11:27:12Z  Alun
pre_adjust() methods added
Revision 1.12  2008/10/22 16:42:02Z  Alun
Early November checkin
Completely reworked pool and Difference_Tracker
Ported to Darwin
Revision 1.11  2008/09/27 17:55:14Z  Alun
Final version built using Indexer.
Revision 1.6  2007/11/15 22:58:12Z  Alun
*/

#pragma once
#ifndef MAF_WE_INCLUDED
#define MAF_WE_INCLUDED 1
#ifndef MAF_EW_INCLUDED
#include "maf_ew.h"
#endif
#ifndef EQUATION_INCLUDED
#include "equation.h"
#endif
#ifndef MAF_TRACE_INCLUDED
#include "maf_trace.h"
#endif


const unsigned PD_CHECK = 1;
const unsigned PD_GM_DIFFERENCE = 2;
const unsigned PD_DR_DIFFERENCE = 4;
const unsigned PD_BOTH = PD_GM_DIFFERENCE | PD_DR_DIFFERENCE;

/* flags to pass to various functions relating to word-differences.
   (GWD= "grow word-differences")
   These don't really belong here but are here because MAF's automata
   growing code does not know about Difference_Tracker, and
   Difference Tracker does not know about Rewriter_Machine, but both
   know about this header */

const unsigned GWD_PRIMARY_ONLY = 1;
const unsigned GWD_KNOWN_TRANSITIONS = 2;

enum Derivation_Step_Type
{
  DST_Swap,
  DST_Cancellation,
  DST_Left_Balance,
  DST_Right_Balance,
  DST_Left_Reduction,
  DST_Right_Reduction,
  DST_Unbalance,
  DST_Transfer,
  DST_Left_Free_Reduction,
  DST_Right_Free_Reduction
};

class Working_Equation; // Forward declaration required by some compilers

/* What a proof trace needs to know about one operation on an equation.
   pid and pos describe a rewrite, g a balancing generator, leading and
   trailing a cancellation */
struct Trace_Op
{
  unsigned short type; // a Derivation_Step_Type
  int pid;
  int pos;
  Ordinal g;
  Word_Length leading;
  Word_Length trailing;
};

class Full_Derivation : public Derivation
{
  friend class Working_Equation;
  friend class Proof_Trace;
  private:
    struct Derivation_Step
    {
      Derivation_Step *next;
      unsigned short type;
      Word_Length offset;
      union
      {
        Node_ID e;
        struct
        {
          Word_Length leading;
          Word_Length trailing;
        } cancellation;
        Ordinal g;
      } data;
      Derivation_Step(Derivation_Step_Type type_) :
        type(type_),
        next(0)
      {}
    };
    Ordinal_Word original_lhs;
    Ordinal_Word original_rhs;
    Derivation_Step * head;
    Derivation_Step * tail;
    bool logged;
    bool started;
    /* Proof tracing: when traced is set every operation is also recorded
       in trace_ops, and trace_base is the derivation of original_rhs from
       original_lhs; see maf_trace.cpp */
    bool traced;
    std::vector<Trace_Op> trace_ops;
    Trace_Chain trace_base;
    std::string trace_origin;
  public:
    Full_Derivation(const Derivation & derivation,const Alphabet & alphabet) :
      Derivation(derivation),
      head(0),
      tail(0),
      logged(false),
      started(false),
      traced(false),
      original_lhs(alphabet),
      original_rhs(alphabet)
    {}
    ~Full_Derivation()
    {
      clean();
    }
    void reset(const Derivation & derivation)
    {
      clean();
      started = false;
      traced = false;
      trace_ops.clear();
      * (Derivation *) this = derivation;
    }

    void add_cancellation(Word_Length leading,Word_Length trailing)
    {
      if (logged)
      {
        Derivation_Step * step = new Derivation_Step(DST_Cancellation);
        step->data.cancellation.leading = leading;
        step->data.cancellation.trailing = trailing;
        add(step);
      }
      if (traced)
        trace(DST_Cancellation,0,0,0,leading,trailing);
    }
    void add_left_balance(Ordinal g)
    {
      if (logged)
      {
        Derivation_Step * step = new Derivation_Step(DST_Left_Balance);
        step->data.g = g;
        add(step);
      }
      if (traced)
        trace(DST_Left_Balance,0,0,g);
    }
    void add_left_free_reduction(Word_Length offset,Ordinal g)
    {
      if (logged)
        add_free_reduction(DST_Left_Free_Reduction,offset,g);
      if (traced)
        trace(DST_Left_Free_Reduction,0,offset,g);
    }
    void add_left_reduction(Node_List & nl,const Node_Manager &nm)
    {
      if (logged)
        add_reduction(DST_Left_Reduction,nl,nm);
    }
    void add_right_balance(Ordinal g)
    {
      if (logged)
      {
        Derivation_Step * step = new Derivation_Step(DST_Right_Balance);
        step->data.g = g;
        add(step);
      }
      if (traced)
        trace(DST_Right_Balance,0,0,g);
    }
    void add_right_free_reduction(Word_Length offset,Ordinal g)
    {
      if (logged)
        add_free_reduction(DST_Right_Free_Reduction,offset,g);
      if (traced)
        trace(DST_Right_Free_Reduction,0,offset,g);
    }
    void add_right_reduction(Node_List & nl,const Node_Manager &nm)
    {
      if (logged)
        add_reduction(DST_Right_Reduction,nl,nm);
    }
    void add_swap()
    {
      if (logged)
      {
        Derivation_Step * step = new Derivation_Step(DST_Swap);
        add(step);
      }
      if (traced)
        trace(DST_Swap);
    }
    void add_transfer(Word_Length l)
    {
      if (logged)
      {
        Derivation_Step * step = new Derivation_Step(DST_Transfer);
        step->offset = l;
        add(step);
      }
      if (traced)
        trace(DST_Transfer,0,l);
    }
    void add_unbalance(Ordinal g)
    {
      if (logged)
      {
        Derivation_Step * step = new Derivation_Step(DST_Unbalance);
        step->data.g = g;
        add(step);
      }
      if (traced)
        trace(DST_Unbalance,0,0,g);
    }
    void log()
    {
      logged = true;
    }
    void print(Container & container,Output_Stream * stream,Rewriter_Machine &rm);
    void set_original(const Working_Equation &we);
  private:
    void trace(Derivation_Step_Type type,int pid = 0,int pos = 0,
               Ordinal g = 0,Word_Length leading = 0,Word_Length trailing = 0)
    {
      Trace_Op op;
      op.type = (unsigned short) type;
      op.pid = pid;
      op.pos = pos;
      op.g = g;
      op.leading = leading;
      op.trailing = trailing;
      trace_ops.push_back(op);
    }
    void add(Derivation_Step * step)
    {
      if (tail)
        tail->next = step;
      else
        head = step;
      tail = step;
    }
    void add_reduction(Derivation_Step_Type type,Node_List & nl,const Node_Manager &nm);
    void add_free_reduction(Derivation_Step_Type type,Word_Length offset,Ordinal g)
    {
      Derivation_Step * step = new Derivation_Step(type);
      step->offset = offset;
      step->data.g = g;
      add(step);
    }
    void clean();
};

struct Equation_Properties
{
  Ordinal lvalue;
  Ordinal rvalue;
  unsigned short eq_flags;
  Word_Length l;
  Word_Length r;
  Total_Length total;
  Node_Reference e;
  Node_Reference old_e;
  Node_Reference prefix;
  Node_Reference rhs_node;
  Node_Reference trailing_subword;
};

class Working_Equation : public Equation
{
  friend class Rewriter_Machine;
  friend class Proof_Trace;
  BLOCKED(Working_Equation)
  private:
    Node_Manager & nm;
    Equation_Word first_word;
    Equation_Word second_word;
    Equation_Word third_word;
    Equation_Word *lhs;
    Equation_Word *rhs;
    Equation_Word *spare;
    unsigned short eq_flags;
  public:
    bool changed;
    bool swapped;
    bool failed; /* normalise() usually sets this false. It is true afterwards
                    if one or other of the words in the equation could not
                    be reduced without the maximum word length being exceeded */
    char balanced; /* 1 if the equation has been balanced. In this case
                     the original conjunction should be reconsidered.
                      2 if the equation was deliberately unbalanced. */
    Full_Derivation derivation;

    Working_Equation(Node_Manager & nm,
                     const Word & lhs_word);
    Working_Equation(Node_Manager & nm,
                     Equation_Word & lhs_word,
                     Equation_Word & rhs_word,
                     const Derivation &derivation_);
    Working_Equation(Node_Manager & nm,
                     const Word & lhs_word,
                     const Word & rhs_word,
                     const Derivation &derivation_);
    Working_Equation(Node_Manager & nm,Node_Handle lhs,Node_Handle rhs,
                     const Derivation &derivation_,Ordinal transition=-1);
    Working_Equation(Node_Manager & nm,const Derivation &derivation_);
    // adjust_for_dm() is destructive. The equation is not valid after this
    Node_Reference adjust_for_dm(Node_Reference * hstate,Ordinal * multiplier,bool check_balance);
    Node_Reference pre_adjust(Node_Reference * hstate) const
    {
      return pre_adjust(hstate,false);
    }

    Equation_Word & lhs_word()
    {
      return *lhs;
    }
    /* learn() requests that the specified equation be "learned" by the
       Rewriter_Machine. The return value is 1 for success, 2 if the equation
       has been pooled, and 0 if the quation is either unnecessary or has
       been rejected */
    int learn(unsigned flags);
    unsigned learn_differences(Node_Handle e,unsigned flags);
    void log();
    bool normalise(unsigned flags = AE_REDUCE|AE_SAFE_REDUCE);
    bool prepare_for_insert(unsigned flags = 0);
    void print_derivation(Output_Stream * stream,ID pool_id = -1);
    int rhs_changed()
    {
      *spare = *lhs;
      bool failed = false;
      spare->reduce(false,&failed);
      if (failed)
        return -2;
      return rhs->compare(*spare);
    }
    Equation_Word & rhs_word()
    {
      return *rhs;
    }
    bool set_lhs(bool cancel_equal = false)
    {
      log();
      swapped = false;
      int cmp = lhs->compare(*rhs);
      if (cmp <= 0)
      {
        if (!cmp)
        {
          if (cancel_equal)
          {
            lhs->set_length(0);
            rhs->set_length(0);
          }
          return false;
        }
        Equation_Word * temp = lhs;
        lhs = rhs;
        rhs = temp;
        derivation.add_swap();
        swapped = true;
      }
      return true;
    }
    Equation_Word & spare_word()
    {
      return *spare;
    }
    // const methods from now on
    void fatal_error(String message) const;
    virtual const Word & get_lhs() const
    {
      return (const Word &) *lhs;
    }
    virtual const Word & get_rhs() const
    {
      return (const Word &) *rhs;
    }
    const Equation_Word & lhs_word() const
    {
      return *lhs;
    }
    const Equation_Word & rhs_word() const
    {
      return *rhs;
    }
    Total_Length total_length() const
    {
      return lhs->word_length + rhs->word_length;
    }
    Total_Length total_g_length() const;
    void read_properties(Equation_Properties * ep);
  private:
    unsigned cancel_common();
    // const below is a lie when adjust is true
    Node_Reference pre_adjust(Node_Reference * hstate,bool adjust) const;
    bool geodesic_balance(unsigned flags);
    bool balance(unsigned flags);
};

#endif
