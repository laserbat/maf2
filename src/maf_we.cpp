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


// $Log: maf_we.cpp $
// Revision 1.20  2010/06/10 13:57:43Z  Alun
// All tabs removed again
// Revision 1.19  2010/05/29 08:38:49Z  Alun
// Lots of source code moved out to maf_ew.cpp
// Changes required by new Node_Manager interface
// Balancing code reorganised with balancing for geodesic
// and non-geodesic word-orderings now handled by separate methods
// Revision 1.18  2009/11/10 08:15:02Z  Alun
// Changed calculation of initial g_word for coset word-differences.
// Various minor changes to allow for future increase to maximum word length
// Revision 1.17  2009/10/10 09:47:10Z  Alun
// Compatibility of FSA utilities with KBMAG improved.
// Compatibility programs produce fewer FSAs KBMAG would not (if any)
// Revision 1.16  2009/09/16 07:50:29Z  Alun
// Additional source code changes needed to get clean compilation on latest GNU
//
// Revision 1.15  2009/09/13 19:14:55Z  Alun
// printf() style functions cleaned up considerably
// code now compiles cleanly in g++ both 32-bit and 64-bit
// all tabs removed from source files
// Revision 1.14  2008/12/29 20:06:14Z  Alun
// Changes for support of normal subgroup coset systems
// Revision 1.14  2008/11/02 18:57:14Z  Alun
// Early November checkin
// Completely reworked pool and Difference_Tracker
// Ported to Darwin
// Revision 1.13  2008/10/10 08:53:45Z  Alun
// diagnostics removed
// Revision 1.12  2008/09/28 20:48:42Z  Alun
// Final version built using Indexer.
// Revision 1.7  2007/12/20 23:25:42Z  Alun
//
#include "maf_rm.h"
#include "mafnode.h"
#include "maf_nm.h"
#include "container.h"
#include "maf_wdb.h"
#include "maf_we.h"

void Full_Derivation::add_reduction(Derivation_Step_Type type,Node_List & nl,
                                    const Node_Manager &nm)
{
  /* Transfer the reduction history of an Equation_Word into a derivation.
     It is necessary to copy reduction history like this because during
     balancing we may frequently perform a trial reduction that we
     then do not use.
  */

  Node_Reference e;
  int offset;
  while (nl.use(&e,nm,&offset))
  {
    Derivation_Step * step = new Derivation_Step(type);
    step->offset = Word_Length(offset);
    step->data.e = e;
    add(step);
  }
}

/**/

void Full_Derivation::print(Container &container,Output_Stream * stream,Rewriter_Machine & rm)
{
  ID pid1 = 0;
  ID pid2 = 0;
  const Alphabet & alphabet = rm.alphabet();
  const Node_Manager & nm = (const Rewriter_Machine &) rm;
  Ordinal_Word ow1(alphabet);
  Ordinal_Word ow2(alphabet);
  String_Buffer sb1;
  String_Buffer sb2;
  Ordinal g1 = PADDING_SYMBOL,g2 = PADDING_SYMBOL;

  switch (type)
  {
    case BDT_Relator_Multiply:
      Node::fast_find_node(nm,data.overlap.e1)->read_word(&ow1,nm);
      pid2 = rm.get_derivation_id(data.overlap.e2,false);
      break;
    case BDT_Equal_RHS:
    case BDT_Conjunction:
      pid2 = rm.get_derivation_id(data.overlap.e2,false);
      /* fall through */
    case BDT_Right_Conjugate:
    case BDT_Known:
    case BDT_Inversion_Closure:
    case BDT_Double:
      pid1 = rm.get_derivation_id(data.overlap.e1,false);
      break;
    case BDT_Left_Conjugate:
    case BDT_Normal_Closure:
      pid1 = rm.get_derivation_id(data.conjugate.e1,false);
      break;
    case BDT_Pool:
    case BDT_Pool_Conjugate:
    case BDT_Equal_LHS:
      pid1 = data.pool_id+1;
      break;
    case BDT_Reverse:
    case BDT_Update:
    case BDT_Half_Difference_1:
    case BDT_Half_Difference_2:
      if (data.computation.state)
        Node::fast_find_node(nm,data.computation.state)->read_word(&ow1,nm);
      else
        ow1 = original_lhs;
      if (data.computation.old_state)
        Node::fast_find_node(nm,data.computation.old_state)->read_word(&ow2,nm);
      if (data.computation.ti >= 0)
        alphabet.product_generators(&g1,&g2,data.computation.ti);
      break;
    case BDT_Unspecified:
    case BDT_Axiom:
    case BDT_Diff_Reduction:
    case BDT_Trivial:
      break;
  }
  Derivation_Step * step = head;
  bool need_original = step && (step->type != DST_Swap || step->next);
  bool long_format = need_original;

  if (step && type != BDT_Conjunction)
    long_format = true;

  if (long_format)
  {
    container.output(stream,"(\n  ");
    if (need_original)
      container.output(stream,"%s=%s\n  ",original_lhs.format(&sb1).string(),
                                          original_rhs.format(&sb2).string());
  }
  else
    container.output(stream,"  (");

  switch (type)
  {
    case BDT_Axiom:
      container.output(stream,"Axiom");
      break;
    case BDT_Trivial:
      container.output(stream,"Trivially");
      break;
    case BDT_Known:
      container.output(stream,"from equation " FMT_ID ,pid1);
      break;
    case BDT_Conjunction:
      if (step && step->type == DST_Swap && !need_original)
      {
        container.output(stream,"Swapped ");
        Derivation_Step * next = step->next;
        delete step;
        step = next;
      }
      container.output(stream,"conjunction of " FMT_ID "," FMT_ID " at %d",pid1,pid2,
                       data.overlap.offset);
      break;
    case BDT_Equal_RHS:
      if (pid1 == pid2)
        container.output(stream,"from old and new RHS of " FMT_ID ,pid1);
      else
        container.output(stream,"from equal RHS of " FMT_ID "," FMT_ID ,pid1,pid2);
      break;
    case BDT_Right_Conjugate:
      container.output(stream,"conjugate of " FMT_ID,pid1);
      break;
    case BDT_Left_Conjugate:
      container.output(stream,"via left multiplication of " FMT_ID " by %s",pid1,
                       rm.alphabet().glyph(data.conjugate.g).string());
      break;
    case BDT_Normal_Closure:
      container.output(stream,"from normal closure on " FMT_ID " by %s",pid1,
                       rm.alphabet().glyph(data.conjugate.g).string());
      break;

    case BDT_Diff_Reduction:
      container.output(stream,"proved by word-differences ");
      break;
    case BDT_Update:
      if (data.computation.ti == -1)
        container.output(stream,"from new and old values of w^-1 for w=%s",
                         ow1.format(&sb1).string());
      else if (!data.computation.old_state)
        container.output(stream,"(old and new values of x^-1*w*y for w=%s (x,y)=(%s,%s)",
                         ow1.format(&sb1).string(),alphabet.glyph(g1).string(),
                         alphabet.glyph(g2).string());
      else
        container.output(stream,
                         "from w1=w2=>x^-1*w1*y=x^-1*w2*y and w1=%s,w2=%s (x,y)=(%s,%s)",
                         ow1.format(&sb1).string(),ow2.format(&sb2).string(),
                         alphabet.glyph(g1).string(),alphabet.glyph(g2).string());
      break;

    case BDT_Reverse:
      if (data.computation.ti == -1)
        container.output(stream,"(w^-1)^-1=w for w=%s",ow1.format(&sb1).string());
      else
        container.output(stream,"x*(x^-1*w*y)*y^-1=w for w=%s (x,y)=(%s,%s)",
                         ow1.format(&sb1).string(),alphabet.glyph(g1).string(),
                         alphabet.glyph(g2).string());
      break;

    case BDT_Half_Difference_1:
      container.output(stream,
                       "(x^-1*w*y)^-1*x^-1=(w*y)^-1= for w=%s (x,y)=(%s,%s)",
                       ow1.format(&sb1).string(),alphabet.glyph(g1).string(),
                       alphabet.glyph(g2).string());
      break;
    case BDT_Half_Difference_2:
      container.output(stream,
                       "w^-1*x=((x^-1*w*y)*y^-1)^-1 for w=%s (x,y)=(%s,%s)",
                       ow1.format(&sb1).string(),alphabet.glyph(g1).string(),
                       alphabet.glyph(g2).string());
      break;
    case BDT_Inversion_Closure:
      container.output(stream,
                       "via word-difference %d of equation " FMT_ID " and its conjugate",
                       data.overlap.offset,pid1);
      break;
    case BDT_Double:
      container.output(stream,
                       "word-difference %d of equation " FMT_ID " from left and right",
                       data.overlap.offset,pid1);
      break;
    case BDT_Pool:
      container.output(stream,"Old P" FMT_ID ,pid1);
      break;
    case BDT_Pool_Conjugate:
      container.output(stream,"Conjugate of P" FMT_ID ,pid1);
      break;
    case BDT_Equal_LHS:
      container.output(stream,"Old and new RHS for P" FMT_ID ,pid1);
      break;
    case BDT_Relator_Multiply:
      container.output(stream,
                       "u = v -> x*u*v^-1 = x for x=%s and equation " FMT_ID ")",
                       ow1.format(&sb1).string(),pid2);
      break;
    case BDT_Unspecified:
    default:
      MAF_INTERNAL_ERROR(container,
                         ("Unexpected value %d for type in"
                         " Full_Derivation::print()\n",type));
      break;
  }
  if (step)
    container.output(stream,"\n");
  Derivation_Step * next;
  for (;step;step = next)
  {
    next = step->next;
    switch (step->type)
    {
      case DST_Swap:
        container.output(stream,"  Swap LHS and RHS\n");
        break;
      case DST_Cancellation:
        if (step->data.cancellation.leading)
          if (step->data.cancellation.trailing)
            container.output(stream,"  Cancel %d symbols on left and %d on right\n",
                             step->data.cancellation.leading,
                             step->data.cancellation.trailing);
          else
            container.output(stream,"  Cancel %d symbols on left\n",
                             step->data.cancellation.leading);
        else
          container.output(stream,"  Cancel %d symbols on right\n",
                           step->data.cancellation.trailing);
        break;
      case DST_Left_Balance:
        for (;;)
        {
          ow1.set_length(0);
          ow1.append(step->data.g);
          ow2 = ow1 + ow2;
          if (!next || next->type != step->type)
            break;
          delete step;
          step = next;
          next = step->next;
        }
        container.output(stream,"  Left multiply by %s and reduce LHS\n",
                         ow2.format(&sb1).string());
        break;
      case DST_Right_Balance:
      case DST_Unbalance:
        ow1.set_length(0);
        for (;;)
        {
          ow1.append(step->data.g);
          if (!next || next->type != step->type)
            break;
          delete step;
          step = next;
          next = step->next;
        }
        if (step->type == DST_Unbalance)
          container.output(stream,"  Right multiply by %s and reduce RHS\n",
                           ow1.format(&sb1).string());
        else
          container.output(stream,"  Right multiply by %s and reduce LHS\n",
                           ow1.format(&sb1).string());
        break;
      case DST_Left_Reduction:
      case DST_Right_Reduction:
        {
          if (next && next->type==step->type)
            container.output(stream,step->type == DST_Left_Reduction ?
                             "  Apply LHS reductions:" :
                             "  Apply RHS reductions:");
          else
            container.output(stream,step->type == DST_Left_Reduction ?
                             "  Apply LHS reduction:" :
                             "  Apply RHS reduction:");

          bool started = false;
          for (;;)
          {
            if (started)
               container.output(stream,",");
            ID id = rm.get_derivation_id(step->data.e,false);
            container.output(stream," " FMT_ID " at %d",id,step->offset);
            if (!next || next->type != step->type)
              break;
            started = true;
            delete step;
            step = next;
            next = step->next;
          }
          container.output(stream,"\n");
        }
        break;
      case DST_Right_Free_Reduction:
        container.output(stream,"  Apply %s*%s=IdWord at %d on RHS\n",
                         alphabet.glyph(rm.inverse(step->data.g)).string(),
                         alphabet.glyph(step->data.g).string(),step->offset);
        break;
    }
    delete step;
  }
  container.output(stream,")\n");
  head = tail = 0;
}

void Full_Derivation::set_original(const Working_Equation & we)
{
  if (logged && !started)
  {
    original_lhs = we.lhs_word();
    if (type != BDT_Trivial)
      original_rhs = we.rhs_word();
    else
      original_rhs = we.lhs_word();
    started = true;
  }
}

/**/

void Full_Derivation::clean()
{
  Derivation_Step * next;
  for (Derivation_Step * step = head;step;step = next)
  {
    next = step->next;
    delete step;
  }
  head = tail = 0;
}

/**/

void Working_Equation::print_derivation(Output_Stream * stream,ID pool_id)
{
  /* Print out the derivation of an equation. We need to assign ids to
     the equations because the e->id() value is not set up yet, and even
     worse may change */

  Container & container = nm.maf.container;
  if (derivation.type != BDT_Pool || derivation.head || pool_id != derivation.data.pool_id)
  {
    bool force_derivation = false;
    if (pool_id >= 0)
      container.output(stream,"P" FMT_ID ":",pool_id+1);
    else
    {
      Node_ID e = lhs->state[lhs->word_length];
      State_ID id = nm.rm.get_derivation_id(e,true);
      container.output(stream,FMT_ID ":",id);
    }
    if (derivation.type != BDT_Pool || derivation.head || pool_id < 0)
      print(container,stream);
    else
      force_derivation = true;
    if (nm.maf.options.log_flags & LOG_DERIVATIONS || force_derivation)
      derivation.print(container,stream,nm.rm);
  }
}

/**/

Working_Equation::Working_Equation(Node_Manager & nm_,
                                   const Word & lhs_word,
                                   const Word & rhs_word,
                                   const Derivation &derivation_) :
  nm(nm_),
  first_word(nm_,lhs_word),
  second_word(nm_,rhs_word),
  third_word(nm_,lhs_word.length()),
  lhs(&first_word),
  rhs(&second_word),
  spare(&third_word),
  derivation(derivation_,nm_.alphabet()),
  changed(false),
  balanced(0),
  eq_flags(0)
{
}

/**/

Working_Equation::Working_Equation(Node_Manager & nm_,
                                   const Derivation &derivation_) :
  nm(nm_),
  first_word(nm_,1),
  second_word(nm_,1),
  third_word(nm_,1),
  lhs(&first_word),
  rhs(&second_word),
  spare(&third_word),
  derivation(derivation_,nm_.alphabet()),
  changed(false),
  balanced(0),
  eq_flags(0)
{
}

/**/

Working_Equation::Working_Equation(Node_Manager & nm_,
                                   Equation_Word & lhs_word,
                                   Equation_Word & rhs_word,
                                   const Derivation &derivation_) :
  nm(nm_),
  first_word(lhs_word),
  second_word(rhs_word),
  third_word(nm_,lhs_word.word_length),
  lhs(&first_word),
  rhs(&second_word),
  spare(&third_word),
  derivation(derivation_,nm_.alphabet()),
  changed(false),
  balanced(0),
  eq_flags(0)
{
}

/**/

Working_Equation::Working_Equation(Node_Manager & nm_,
                                   const Word & lhs_word) :
  nm(nm_),
  first_word(nm_,lhs_word),
  second_word(nm_,lhs_word),
  third_word(nm_,first_word.word_length),
  lhs(&first_word),
  rhs(&second_word),
  spare(&third_word),
  derivation(Derivation(BDT_Trivial),nm_.alphabet()),
  changed(false),
  balanced(0),
  eq_flags(0)
{
  if (nm.maf.options.log_flags & LOG_DERIVATIONS)
    derivation.log();
}

/**/

Working_Equation::Working_Equation(Node_Manager & nm_,
                                   Node_Handle lhs_state,
                                   Node_Handle rhs_state,
                                   const Derivation &derivation_,
                                   Ordinal transition) :
  nm(nm_),
  first_word(nm_,lhs_state,transition),
  second_word(nm_,rhs_state),
  third_word(nm_,lhs_state->length(nm_)),
  lhs(&first_word),
  rhs(&second_word),
  spare(&third_word),
  derivation(derivation_,nm_.alphabet()),
  changed(false),
  balanced(0),
  eq_flags(0)
{
}

/**/

Working_Equation::Working_Equation(const Working_Equation & other) :
  Equation(other), // needed to shut up g++ - does nothing
  nm(other.nm),
  first_word(*other.lhs),
  second_word(*other.rhs),
  third_word(other.nm,other.lhs->word_length),
  lhs(&first_word),
  rhs(&second_word),
  spare(&third_word),
  derivation(other.derivation,other.nm.alphabet()),
  changed(false),
  balanced(0),
  eq_flags(0)
{
}

/**/

unsigned Working_Equation::cancel_common()
{
  /* Removes common cancellable symbols from start and end of LHS and RHS
     return value is 1 if there is cancellable suffix
                     2 if there is cancellable prefix
                     3 if both are present
  */

  Rewriter_Machine &rm = nm.rm;
  Word_Length l = lhs->word_length;
  Word_Length r = rhs->word_length;
  unsigned retcode = 0;

  Word_Length start = 0;
  Word_Length common = r <= l ? r : l;
  Word_Length trailing = 0;
  /* see if there is a common cancellable suffix */
  {
    Ordinal rvalue;
    Word_Length i;
    for (i = 1; i <= common;i++)
    {
      if ((rvalue = lhs->values[l-i]) != rhs->values[r-i] ||
          !rm.right_cancels(rvalue))
        break;
    }
    if (i != 1)
    {
      trailing = i-1;
      l -= trailing;
      r -= trailing;
      common -= trailing;
      retcode = 1;
    }
  }
  /* see if there is a common cancellable prefix */
  {
    Ordinal lvalue;
    while (start < common && (lvalue = lhs->values[start]) == rhs->values[start] &&
           rm.left_cancels(lvalue))
      start++;
    if (start)
      retcode |= 2;
  }
  if (retcode)
  {
    /* Complete cancellations */
    lhs->slice(start,l);
    rhs->slice(start,r);
    derivation.add_cancellation(start,trailing);
    changed = true;
  }
  return retcode;
}

/**/

static Total_Length cancellable_length(const Rewriter_Machine & rm,const Word & lhs,const Word & rhs)
{
  /* Counts the common cancellable symbols from start and end of LHS and RHS */
  Word_Length l = lhs.length();
  Word_Length r = rhs.length();
  const Ordinal * lvalues = lhs.buffer();
  const Ordinal * rvalues = rhs.buffer();
  Total_Length retcode = 0;

  Word_Length start = 0;
  Word_Length common = r <= l ? r : l;
  /* see if there is a common cancellable suffix */
  Word_Length i;
  for (i = 1; i <= common;i++)
    if (lvalues[l-i] != rvalues[r-i] || !rm.right_cancels(lvalues[l-i]))
      break;
  if (--i)
  {
    common -= i;
    retcode = i*2;
  }
  /* see if there is a common cancellable prefix */
  while (start < common && lvalues[start] == rvalues[start] &&
         rm.left_cancels(lvalues[start]))
    start++;
  retcode += start*2;
  return retcode;
}

/**/

void Working_Equation::log()
{
  if (nm.maf.options.log_flags & LOG_DERIVATIONS)
  {
    derivation.log();
    derivation.set_original(*this);
  }
}

/**/

bool Working_Equation::normalise(unsigned flags)
{
  /* On entry lhs and rhs are known or stated to be equal.
     We we perform any permitted cancellations and, if necessary,
     transfer elements from lhs to rhs (the lengths should be as equal
     as possible subject to lexical considerations that would interchange
     lhs and rhs).

     The return value is false if the equation is trivial, or to be
     rejected on other grounds, or true if it looks OK.

     When words are not ordered by shortlex it is possible for normalisation
     to fail because one or other of the words in the equation cannot be
     reduced without exceeding the maximum permitted word length. In this
     case MAF will terminate with an error unless the AE_SAFE_REDUCE flag
     flag is passed, whrn it will instead set the "failed" flag in the equation
     and return true, (so that add_equation() can put the equation in the pool
     rather than discarding it).
  */

  eq_flags = derivation.type == BDT_Axiom ? EQ_AXIOM : 0;
  failed = false;
  bool *no_fail = flags & AE_SAFE_REDUCE ? &failed : 0;
  if (nm.pd.is_short)
    no_fail = 0;

  log();
  if (flags & AE_PRE_CANCEL)
    cancel_common();

  if (flags & AE_REDUCE)
  {
    if (rhs->reduce(0,no_fail))
    {
      if (!(flags & AE_COPY_LHS))
        changed = true;
      if (flags & AE_SECONDARY && rhs->word_length &&
          rhs->values[0] == lhs->values[0])
        return false; /* In this case we are checking partial reductions and
                         can ignore equations which have cancellations */
      derivation.add_right_reduction(rhs->reductions,nm);
    }
    else if (flags & AE_SECONDARY)
      return false;
    if (lhs->reduce(flags,no_fail))
    {
      if (flags & AE_SECONDARY)
        return false;
      changed = true;
      derivation.add_left_reduction(lhs->reductions,nm);
    }
  }

  /* ensure lhs is later than rhs */
  if (!set_lhs(true))
    return false;
  if (swapped || changed)
    flags &= ~AE_RIGHT_CONJUGATE;
  if (nm.pd.presentation_type == PT_Monoid_Without_Cancellation)
    return true;
  if (flags & AE_NO_BALANCE)
  {
    cancel_common();
    return true;
  }
  return nm.pd.is_short ? geodesic_balance(flags) : balance(flags);
}

/**/

bool Working_Equation::geodesic_balance(unsigned flags)
{
  /* Since we are only called from normalise(), upon entry we are sure
     the lhs and rhs are different and the lhs is different.
     We try to reduce the length of the lhs by multiplying on the left or
     right by the inverse of the appropriate generators, while
     keeping lhs > rhs */
  Rewriter_Machine &rm = nm.rm;

  for (;;)
  {
    if (cancel_common())
      flags &= ~(AE_KEEP_LHS+AE_SECONDARY+AE_RIGHT_CONJUGATE);
    Word_Length l = lhs->word_length;
    Word_Length r = rhs->word_length;

    if (rhs->language != language_L0)
    {
      /* doesn't mean rhs is not reduced, just we don't know whether it is */
      if (rhs->reduce())
      {
        derivation.add_right_reduction(rhs->reductions,nm);
        changed = true;
        continue;
      }
      r = rhs->word_length;
    }

    /* Now we try to equalise left and right hand sides by removing
       letters from end of lhs and replacing them with their inverses on
       the end of the rhs (i.e. right multiplication by inverse of
       rightmost element of LHS).
       We should not try making the RHS longer than the LHS on the off
       chance it will reduce to a lesser word, because this is something
       it is appropriate to discover when right conjugating.
    */

    Ordinal inverse_value;

    while (l > r + 2 &&
           (inverse_value = rm.inverse(lhs->values[l-1])) != INVALID_SYMBOL)
    {
      lhs->set_length(--l);
      rhs->append(inverse_value);
      r++;
      changed = true;
      flags &= ~AE_KEEP_LHS;
      derivation.add_right_balance(inverse_value);
    }

    if (l == r+2 && r &&
        (inverse_value = rm.inverse(lhs->values[l-1])) != INVALID_SYMBOL &&
        (!nm.pd.is_shortlex || lhs->values[0] > rhs->values[0]))
    {
      rhs->append(inverse_value);
      if (nm.pd.is_shortlex || Subword(*lhs,0,l-1).compare(*rhs) > 0)
      {
        lhs->set_length(--l);
        r++;
        changed = true;
        flags &= ~AE_KEEP_LHS;
        derivation.add_right_balance(inverse_value);
      }
      else
        rhs->set_length(r);
    }

    if (!(flags & AE_KEEP_LHS) && l == r+2 &&
        (inverse_value = rm.inverse(lhs->values[0])) != INVALID_SYMBOL)
    {
      /* If lengths still differ by 2 we may be able to move an element
         from the left of the lhs without interchanging lhs and rhs.
         Always to do this causes problems for the word-difference
         machines.
         Suppose we have arrived at an equation cdef=da.
         Left balancing changes this to def=Cda.
         As a result we never form a conjugate of the equation with an
         LHS beginning with cd, but the equation cdef=da is also needed
         because it has an irreducible prefix. This equation might be
         the only one for which Cd arise as a word-difference.
         On the other hand we should form the balanced equation first,
         since this will give more simplifications.
         The unbalanced equations are not in the minimal reduction
         engine, but having them speeds up reduction. Since we
         record the language each word belongs to we can generate
         a number of different variations for a reduction engine.
      */
      bool rhs_changed = false;
      spare->allocate(r+1,false);
      spare->set_length(0);
      spare->append(inverse_value);
      spare->join(*spare,*rhs);
      if (spare->reduce())
        rhs_changed = true;
      else
      {
        //if we can't reduce a word containing the RHS as a suffix
        // we certainly can't reduce the RHS, so it is reduced.
        rhs->language = language_L0;
      }
      if (Subword(*lhs,1,l).compare(*spare) > 0)
      {
        Equation_Word * temp = rhs;
        rhs = spare;
        spare = temp;
        r = rhs->word_length;
        l = lhs->slice(1,l);
        changed = true;
        derivation.add_left_balance(inverse_value);
        if (rhs_changed)
        {
          derivation.add_right_reduction(rhs->reductions,nm);
          // we might have reduced rhs to something with common prefix
          // or suffix again
          continue;
        }
      }
    }

    if (!(flags & AE_REDUCE) || rhs->language == language_L0 ||
        !rhs->reduce())
    {
      derivation.add_right_reduction(rhs->reductions,nm);
      return true;
    }
    /* record the reduction now, before the cancellations it may allow */
    derivation.add_right_reduction(rhs->reductions,nm);
    changed = true;
  }
}

/**/

bool Working_Equation::balance(unsigned flags)
{
  /* Since we are only called from normalise(), upon entry we are sure
     the lhs and rhs are different and the lhs is different.
     We try to reduce the length of the lhs by multiplying on the left or
     right by the inverse of the appropriate generators, while
     keeping lhs > rhs */
  Rewriter_Machine &rm = nm.rm;
  MAF & maf = nm.maf;
  int allow_transfer = maf.options.balancing_flags;
  Total_Length limit = rm.status().visible_limit;
  /* With recursive type orderings equations can get extremely long,
     and we may be able to perform very many balancing operations.
     In such cases the equation will almost always be thrown away or
     pooled, and there is no point in balancing it. So we use ceiling
     as a threshold on the sides: if both are longer we will stop balancing*/
  Word_Length ceiling = limit/2 + 8;
  if (ceiling < 60)
    ceiling = 60;
  bool spare_failed = false;
  bool *no_fail = flags & AE_SAFE_REDUCE ? &failed : 0;

  for (;;)
  {
    if (cancel_common())
    {
      flags &= ~(AE_KEEP_LHS+AE_SECONDARY+AE_RIGHT_CONJUGATE);
      allow_transfer = maf.options.balancing_flags;
    }
    Word_Length l = lhs->word_length;
    Word_Length r = rhs->word_length;
    bool re_reduced = false;

    if (rhs->language != language_L0)
    {
      /* doesn't mean rhs is not reduced, just we don't know whether it is */
      if (rhs->reduce(0,no_fail))
      {
        derivation.add_right_reduction(rhs->reductions,nm);
        changed = true;
        continue;
      }
      r = rhs->word_length;
    }
    if (failed)
      return true;

    if (l+r > 2*limit && flags & AE_VERY_DISCARDABLE)
      return false;

    if (nm.pd.is_coset_system)
    {
      /* In this part of the code we do two things:
         We check if the equation contains the coset symbol, and if
         so try to move all H-generators from the LHS to the RHS.

         If all the G-generators are at the same level, or the equation is
         an H-equation and all the H-generators are at the same level we
         do right balancing as we would in the shortlex case since this will
         be much faster than the generic non-geodedic balancing.

         We do not check whether the LHS and RHS of equations obey the
         rules for coset systems specified in Presentation::add_axiom()
         because it is impossible for them not to, and the check might take
         a long time in the case where H-generators are present. However,
         since we need to know the position of the coset symbol in coset
         equations we do perform sanity checks in that case.
      */
      Word_Length ls = INVALID_LENGTH,rs = INVALID_LENGTH; // position of coset symbol

      Word_Length lpl = l; // lpl and rpl will represent the part of the
      Word_Length rpl = r; // equation that we can right balance
      Ordinal llvalue,rlvalue;

      if (lhs->values[0] >= nm.pd.coset_symbol &&
          lhs->values[l-1] < nm.pd.coset_symbol)
      {
        if (lhs->values[0] > nm.pd.coset_symbol)
        {
          ls = l-1;
          while (lhs->values[ls] < nm.pd.coset_symbol)
            ls--;
        }
        else
          ls = 0;
        if (lhs->values[ls] != nm.pd.coset_symbol)
          fatal_error("LHS of coset equation must contain coset symbol!\n");
        if (r == 0)
          fatal_error("RHS of coset equation must not be empty!\n");

        if (rhs->values[0] > nm.pd.coset_symbol)
        {
          rs = r-1;
          while (rhs->values[rs] < nm.pd.coset_symbol)
            rs--;
        }
        else
          rs = 0;
        if (rhs->values[rs] != nm.pd.coset_symbol)
          fatal_error("RHS of coset equation must contain coset symbol!\n");
        lpl -= ls + 1;
        rpl -= rs + 1;

        /* if the equation involves the coset symbol, then if possible
           change equations of the form j*_H*u == k*H*v to
                                          _H*u == J*k*H*v.
            In the case where G is a group  u and v must be different
            (otherwise they would have been cancelled), and u > v,
            But if G is only a monoid we might have u=v and in this case
            we must not balance.
        */
        if (nm.pd.presentation_type==PT_Coset_System_With_Inverses && ls &&
            (nm.pd.g_is_group || lpl != rpl ||
             words_differ(lhs->values + ls+1,rhs->values + rs+1,lpl)))
        {
          derivation.add_transfer(ls);
          rm.invert(spare,*lhs,ls);
          if (r + ls > MAX_WORD)
          {
            failed = true;
            return true;
          }
          rhs->join(*spare,*rhs);
          l = lhs->slice(ls,l);
          r = rhs->word_length;
          ls = 0;
          changed = true;
          rs += ls;
        }
        llvalue = lhs->values[ls+1];
        rlvalue = rhs->values[rs+1];
      }
      else
      {
        llvalue = lhs->values[0];
        rlvalue = rhs->values[0];
      }

      if (llvalue < nm.pd.coset_symbol ? nm.pd.g_level : nm.pd.h_level)
      {
        Ordinal inverse_value;
        while ((lpl > rpl + 2 || lpl >= rpl+2 && rpl && llvalue > rlvalue) &&
               (inverse_value = maf.inverse(lhs->values[l-1])) !=
                INVALID_SYMBOL)
        {
          lhs->set_length(--l);
          rhs->append(inverse_value);
          derivation.add_right_balance(inverse_value);
          r++;
          lpl--;
          if (!rpl)
            rlvalue = inverse_value;
          rpl++;
          changed = true;
          allow_transfer = 2; // we have done all the right-balancing we can
        }
      }

      if (rhs->language != language_L0)
      {
        if (rhs->reduce(0,no_fail))
        {
          derivation.add_right_reduction(rhs->reductions,nm);
          continue;
        }
        if (failed)
          return true;
      }

      /* We may as well continue into non-shortlex balancing */
    }

    /* In the non shortlex case balancing is more complicated because
       rhs can be longer than lhs. But if lhs stays >= rhs we should be
       moving letters from end or start of lhs to rhs.
       Although theoretically we should ignore the word length when
       making this comparison, in practice it seems much better not to
       create "balanced" equations that are far from shortlex balanced.
       If they are needed conjugation will eventually force them to be
       created.
       However it is worth trying to move generators from the LHS to
       the RHS even when the RHS is longer, because the RHS might
       then reduce to a more equal word.
    */
    if (l > ceiling && r > ceiling)
      if (flags & AE_DISCARDABLE && l < r)
        allow_transfer = 0;
      else
        allow_transfer &= 1; /* I would like to put allow_transfer=0 here
                                but this is not safe, though it might be
                                if we could be sure the equation will be
                                thrown away. */
    Ordinal inverse_value;
    bool spare_ok = false;
    while (l > 1 && !re_reduced && allow_transfer & 1 &&
           !(flags & AE_RIGHT_CONJUGATE && !swapped && !changed) &&
           (inverse_value = maf.inverse(lhs->values[l-1]))!=INVALID_SYMBOL &&
           r < MAX_WORD)
    {
      if (!spare_ok)
        *spare = *rhs;
      spare->append(inverse_value);
      if (Subword(*lhs,0,l-1).compare(*spare) <= 0)
      {
        allow_transfer &= 2;
        re_reduced = false;
        break;
      }

      if (spare->language != language_L0)
      {
        spare_failed = false;
        if (spare->reduce(0,&spare_failed))
          re_reduced = true;
        if (spare_failed)
        {
          allow_transfer &= 2;
          break;
        }
      }
      else
        spare_ok = true;
      Word_Length sl = spare->word_length;

      if (maf.options.right_balancing & 2)
      {
        /* We check here that we aren't going to stop an equation being
           inserted when it would have been before balancing */
        if (l+r <= limit && l-1+sl > limit)
        {
          Total_Length nl = l-1+sl-cancellable_length(rm,Subword(*lhs,0,l-1),*spare);
          if (nl > limit || nl==0)
          {
            re_reduced = false;
            break;
          }
        }
      }

      /* If the right balancing style mask has 4 set:
         we check the LHS is still longer than the RHS,
         or that the new RHS is no more over-long than the old one,
         or that the equation is still insertable */
      if (!(maf.options.right_balancing & 4) ||
          sl <= l-1 || sl+1 <= r || sl+l-1 <= limit)
      {
        lhs->set_length(--l);
        Equation_Word * temp = rhs;
        rhs = spare;
        spare = temp;
        if (spare_ok)
          spare->append(inverse_value);
        r = sl;
        derivation.add_right_balance(inverse_value);
        if (re_reduced)
        {
          derivation.add_right_reduction(rhs->reductions,nm);
          balanced = 1;
        }
        else if (l > ceiling && r > ceiling &&
                 flags & AE_DISCARDABLE && l < r)
          allow_transfer = 0;
        changed = true;
        flags &= ~AE_KEEP_LHS;
      }
      else
      {
        re_reduced = false;
        break;
      }
    }

    /* Now we balance from the left */
    while (l > 1 && !(flags & AE_KEEP_LHS) && !re_reduced &&
           allow_transfer & 2 &&
           (inverse_value = maf.inverse(lhs->values[0])) != INVALID_SYMBOL &&
           r < MAX_WORD)
    {
      spare->allocate(r+1,false);
      spare->set_length(0);
      spare->append(inverse_value);
      spare->join(*spare,*rhs);
      /* we don't check that we don't need a reduction to have the
         new RHS less than the new LHS. That seems to cause trouble,
         possibly because it sometimes takes a while for equations
         to get conjugated. */

      spare_failed = false;
      if (spare->reduce(0,&spare_failed))
        re_reduced = true;
      if (spare_failed)
      {
        allow_transfer &= 1;
        break;
      }
      Word_Length sl = spare->word_length;
      if (l+r <= limit && l-1+sl > limit && maf.options.left_balancing < 3)
      {
        Total_Length nl = l-1+sl-cancellable_length(rm,Subword(*lhs,1,l),
                                                    *spare);
        if (nl > limit || nl==0)
        {
          re_reduced = false;
          break;
        }
      }
      if ((sl <= l-1 || sl+1 <= r || maf.options.left_balancing >= 2) &&
          Subword(*lhs,1,l).compare(*spare) > 0)
      {
        Equation_Word * temp = rhs;
        rhs = spare;
        spare = temp;
        r = sl;
        l = lhs->slice(1,l);
        derivation.add_left_balance(inverse_value);
        if (re_reduced)
        {
          derivation.add_right_reduction(rhs->reductions,nm);
          balanced = 1;
        }
        else if (l > ceiling && r > ceiling && flags & AE_DISCARDABLE && l < r)
          allow_transfer = 0;
        changed = true;
      }
      else
      {
        re_reduced = false;
        break;
      }
    }

    if (re_reduced)
      continue;

    if (!(flags & AE_REDUCE) || rhs->language == language_L0 ||
        !rhs->reduce(0,no_fail) || failed)
    {
      derivation.add_right_reduction(rhs->reductions,nm);
      return true;
    }
    /* record the reduction now, before the cancellations it may allow */
    derivation.add_right_reduction(rhs->reductions,nm);
    changed = true;
  }
}

/**/

Node_Reference Working_Equation::pre_adjust(Node_Reference * hstate,bool adjust) const
{
  /* This method calculates the appropriate initial state to use
     for the word-differences of the specified equation.
     If the adjust parameter is true, the method being used by adjust_for_dm()
     to alter the equation so that is suitable for the Difference_Tracker
     and realises the initial state.
     If adjust is false then the method does not disturb the equation and
     does not realise the initial state */

  const Presentation_Data &pd = nm.pd;
  Node_Reference initial_state = nm.start();
  *hstate = initial_state;

  if (pd.is_coset_system)
  {
    Ordinal g = lhs->first_letter();
    if (g == pd.coset_symbol)
    {
      // it doesn't matter if total_length() > MAX_WORD below, though
      // we will get an internal error shortly afterwards
      Equation_Word w0(nm,Word_Length(total_length()));
      Word_Length r = rhs->word_length;
      Word_Length i;
      const Ordinal * right_values = rhs->values;
      for (i = 0; i < r;i++)
        if (right_values[i] == pd.coset_symbol)
          break;
      if (i == r)
        fatal_error("Coset symbol missing on RHS of equation\n");
      Word_Length h_length = i;
      if (pd.presentation_type == PT_Simple_Coset_System)
      {
        /* In this case we don't have any information in the equation
           about the hword, but since u=hv h=uv-1.
           We record the hstate as _H, since the tracker needs this
           to be different from the ordinary initial state */
        w0 = Subword(*rhs,0,1);
        *hstate = w0.realise_state();
      }
      else
      {
        w0 = *rhs;
        w0.slice(0,h_length);
        if (adjust)
          *hstate = w0.realise_state();
        else
        {
          w0.reduce();
          *hstate = w0.state_word(true);
        }
        /* Here we used to calculate the g_word from the h word,
           but elsewhere we calculate it from differences, so this
           is not a very good plan, especially given that h words
           tend to be extremely long */
//        if (!w0.g_word(right_values,h_length))
//          return 0;
      }

      w0 = Subword(*rhs,h_length+1);
      nm.rm.invert(&w0,w0);
      w0 = Subword(*lhs,1) + w0;

      if (adjust)
      {
        initial_state = w0.realise_state();
        lhs->slice(1,WHOLE_WORD);
        rhs->slice(h_length+1,r);
      }
      else
        initial_state = w0.state_word(true);
    }
    else if (g > pd.coset_symbol)
      return Node_Reference(0,0); // We don't care about equations in the H generators alone
  }
  return initial_state;
}

/**/

Node_Reference Working_Equation::adjust_for_dm(Node_Reference * hstate,
                                               Ordinal * multiplier,
                                               bool ensure_balanced)
{
  /* adjust_for_dm() makes the changes required to equations so that they
     can be processed for word-differences. For ordinary equations this
     just involves discovering the last symbol on the left hand side and
     ensuring any necessary padding characters are present. The last symbol
     is not replaced with padding: the difference tracker ignores it as
     required, and it is often convenient to have the symbol still there.

     In the coset case it removes coset/subgroup generators from the equation,
     leaving an untrue G equation in its place, and returns the starting
     state that makes it true.
     i.e. _H*ug=h*_H*v is turned into ug=v and the state of the G word equal
     to h is returned. It can then be scanned by something that processes
     word-differences in a simple-minded way. The H word is also returned,
     because coset word reduction needs to know this. (KBMAG seems to
     put a G word on the left of the separator, but this is no good in MAF
     because the Difference_Tracker and the Rewriter_Machine must have the
     same vocabulary.
     The return value is 0 if the equation is in H generators only
  */

  Node_Reference initial_state = pre_adjust(hstate,true);
  if (!initial_state.is_null())
  {
    const Word_Length l = lhs->word_length;
    const Word_Length r = rhs->word_length;
    *multiplier = lhs->last_letter();
    if (!nm.pd.is_short && ensure_balanced && nm.maf.options.right_balancing != 1)
    {
      /* There is no point calculating word-differences for equations that
         have not been properly balanced as they are certain to be
         eliminated */
      rhs->append(nm.inverse(lhs->last_letter()));
      bool unbalanced = Subword(*lhs,0,Word_Length(l-1)).compare(*rhs) > 0;
      rhs->set_length(r);
      if (unbalanced)
        return Node_Reference(0,0);
    }
    Word_Length m = max(l,r);
    lhs->pad(m);
    rhs->pad(m);
  }
  return initial_state;
}

/**/

void Working_Equation::fatal_error(String message) const
{
  print(nm.maf.container,nm.maf.container.get_stderr_stream());
  MAF_INTERNAL_ERROR(nm.maf.container,("%s",message.string()));
}

/**/

bool Working_Equation::prepare_for_insert(unsigned ae_flags)
{
  failed = false;
  if (!normalise(ae_flags|AE_REDUCE|AE_SAFE_REDUCE))
    return false;

  if (!failed)
  {
    /* normalise may return equation with words without state, for
       example if it performed cancellations after reduction */
    lhs->reduce(AE_KEEP_LHS);
    derivation.add_left_reduction(lhs->reductions,nm);
    rhs->reduce(0);
    derivation.add_right_reduction(rhs->reductions,nm);
  }
  return true;
}

/**/

void Working_Equation::read_properties(Equation_Properties * ep)
{
  /* This method won't work as intended unless the equation has been prepared
     with prepare_for_insert() */
  ep->lvalue = lhs->first_letter();
  ep->rvalue = lhs->last_letter();
  ep->eq_flags = eq_flags;
  ep->l = lhs->word_length;
  ep->r = rhs->word_length;
  ep->total = ep->l+ep->r;
  if (!failed)
  {
    ep->old_e = lhs->state_lhs(false);
    if (ep->old_e->length(nm) == lhs->word_length)
      ep->e = ep->old_e;
    else
      ep->e = Node_Reference(0,0);
    ep->trailing_subword = ep->old_e->primary(nm,ep->old_e);
    if (ep->trailing_subword == ep->e)
      if (ep->e->is_final())
        ep->trailing_subword = ep->trailing_subword->overlap_suffix(nm);
      else
        ep->trailing_subword = Node_Reference(0,0);
  }
  else
    ep->old_e = ep->e = ep->trailing_subword = Node_Reference(0,0);
  ep->prefix = lhs->state_prefix();
  ep->rhs_node = rhs->state_word(true);
}

/**/

Total_Length Working_Equation::total_g_length() const
{
  Total_Length answer = total_length();
  if (nm.pd.is_coset_system && lhs->values[0] == nm.pd.coset_symbol)
  {
    Total_Length h_length = 0;
    const Ordinal * rvalues = rhs->values;
    for (h_length = 0; rvalues[h_length] >= nm.pd.coset_symbol;h_length++)
      ;
    h_length++;
    answer -= h_length;
  }
  return answer;
}

/**/

int Working_Equation::learn(unsigned flags)
{
  return nm.rm.add_equation(this,flags);
}
