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


// $Log: maf_rm.cpp $
// Revision 1.25  2011/05/20 10:17:01Z  Alun
// detabbed
// Revision 1.24  2010/07/02 13:57:28Z  Alun
// Behaved badly with finite index subgroups when detect_finite_index was not on
// Revision 1.23  2010/06/18 00:14:28Z  Alun
// Pointless code in coset_extras() removed.
// meaning of various special_overlaps options changed
// Revision 1.22  2010/06/10 13:57:39Z  Alun
// All tabs removed again
// Revision 1.21  2010/05/30 07:56:44Z  Alun
// June 2010 version.
// New style Node_Manager, changes to options to remove "magic number" options
// Changes to handling of coset systems to try to improve stability with
// named subgroup generators, and to consider overlaps that still need to be
// considered after finite index detected
// Revision 1.20  2009/11/11 00:32:19Z  Alun
// normalised_limit renamed to visible_limit
// full conjugation did not work with coset systems
// improved stability of no_pool option
// various new options, including more aggressive version of detect_finite_index
// Revision 1.19  2009/10/13 22:21:23Z  Alun
// Compatibility of FSA utilities with KBMAG improved.
// Compatibility programs produce fewer FSAs KBMAG would not (if any)
// Revision 1.18  2009/09/14 09:44:09Z  Alun
// printf() style functions cleaned up considerably
// code now compiles cleanly in g++ both 32-bit and 64-bit
// all tabs removed from source files
// Warn if word-differences requested for order which does not support automation.
// Revision 1.17  2009/09/02 07:21:38Z  Alun
// Long comment drastically changed.
// explore_dm() must try to ensure normal subgroup coset system is really
// a normal subgroup by conjugating the Schreier generators.
// "aborted" and magic numbers replaced by rm_state and named constants.
// no_pool option was broken
// Revision 1.16  2009/08/23 19:16:27Z  Alun
// Support for normal subgroup coset systems added
// Revision 1.19  2008/11/05 02:39:57Z  Alun
// Economised construction of temporary coset word-acceptor
// Revision 1.18  2008/11/04 23:26:49Z  Alun
// Economised construction of temporary coset word-acceptor
// Revision 1.17  2008/11/03 11:39:40Z  Alun
// Early November checkin
// Completely reworked pool and Difference_Tracker
// Ported to Darwin
// Revision 1.15  2008/10/01 01:37:37Z  Alun
// diagnostics removed
// Revision 1.14  2008/09/29 22:43:48Z  Alun
// Switch to using Hash+Collection_Manager instead of Indexer
// Revision 1.7  2007/12/20 23:25:42Z  Alun
//

/* This module implements MAF's Rewriter_Machine class and associated
   classes, and therefore is, to a large extent, the heart of the entire
   MAF package. A Rewriter_Machine object is very complex, so
   Rewriter_Machine delegates various aspects of its work to several other
   classes, most notably Node_Manager (member nm), which is exposed directly
   to several other parts of MAF, and Job_Manager (jm), which is a nested
   class that is responsible for many scheduling and performing many of the
   updates to the Rewriter_Machine.

   Most of the data that constitutes the rewriting system is stored and
   managed by Node_Manager, so the rest of the very long comment that used
   to be here has been moved to maf_nm.cpp.

   Several other members of Rewriter_Machine are of considerable importance:

     "pool" contains a database of equations that we have discovered but
     chosen not to keep in the tree managed by Node_Manager. The "pool" is
     periodically examined and its contents are used to update Node_Manager.

     "em" is Equation_Manager. All equations are added via this class.

     "dt" is Difference_Tracker, and is responsible for computing word
     differences.

   Rewriter_Machine does three very important things:

   i) It contains MAF's implementation of the Knuth Bendix procedure, which
   is the method by which an initial rewriting system containing axioms
   from some mathematical object, is turned into a rewriting system with
   a proper associative multiplication.

   ii) It is responsible for maintaining all the information in the rewriting
   system. The rewriting system is needed even when we are not going to
   perform the Knuth-Bendix procedure, as its data is also used when
   creating an automatic structure.

   iii) It is responsible for making many strategic and tactical decisions
   which will affect our chances of success.

   Rewriter_Machine's implementation has some special features which might
   be different from other implementations.

   1) Whenever a new equation is inserted into the system, all equations for
   which a prefix of the LHS is now reducible are removed from the tree.
   All such equations are kept in a list which is examined by Job_Manager
   from time to time, since the reduction of the LHS and the old RHS may
   give us more new equations.

   2) A Rewriter_Machine will almost always have some equations in which a
   suffix of the LHS is reducible in the tree. Such equations, which are
   called "secondary equations" are included when the LHS and RHS have no
   cancellable common suffix or prefix. In a "minimal" rewriting system
   equations are only included when no prefix or suffix of the LHS is
   reducible. However in such a machine reduction takes longer, because it
   is necessary to back-track more. Also, at least some of the additional
   equations that Rewriter_Machine creates are usually needed for the
   calculation of the word-differences needed for a general multiplier.

   Equations are always inserted one at a time, and the modifications
   needed to satisfy the minimum consistency requirements are made after each
   insertion.
*/

#include <time.h>
#include <string.h>
#include "fsa.h"
#include "maf_rm.h"
#include "maf_nm.h"
#include "maf_jm.h"
#include "mafnode.h"
#include "equation.h"
#include "keyedfsa.h"
#include "container.h"
#include "variadic.h"
#include "maf_trace.h"
#include "maf_dt.h"
#include "maf_dr.h"
#include "maf_ont.h"
#include "maf_wdb.h"
#include "maf_rws.h"
#include "maf_em.h"

const int Priority_Extended = -100;
const int Priority_Consider = -101;
const int Priority_General = -102;
const int Priority_Critical = 0;
const Total_Length POOL_DEPTH = 8;

const unsigned RMS_Normal = 0;
const unsigned RMS_Aborted_Pass = 1;
const unsigned RMS_Hurrying = 2;
const unsigned RMS_Too_Tight = 4;
const unsigned RMS_Aborting = 8;
const unsigned RMS_Abnormal = 13;
const unsigned RMS_Starting = 16;
const unsigned RMS_Limited = 32;
const unsigned RMS_Widen = 64;
const unsigned RMS_Stick = 128;
const unsigned RMS_Want_Finite_Check = 256;

/**/

void Rewriter_Machine::remove_differences(Equation_Handle e)
{
  if (dt)
  {
    dt->remove_differences(e);
    e->node(nm).clear_flags(EQ_HAS_DIFFERENCES|EQ_HAS_PRIMARY_DIFFERENCES|EQ_INTERESTING);
  }
}

/**/

bool Rewriter_Machine::is_interesting_equation(const Working_Equation & we,unsigned flags)
{
  if (stats.want_differences && dt)
  {
    if (!(flags & AE_CORRECTION+AE_SECONDARY))
    {
      if (!stats.favour_differences)
        return false;
      if (!maf.options.assume_confluent && !pd.g_level)
      {
        if (we.lhs->word_length < we.rhs->word_length ||
            we.total_length() > stats.discard_limit &&
            we.total_length() > nm.overlap_filter.forbidden_limit)
          return false;
        if (pd.inversion_difficult && we.total_length() > stats.visible_limit*2)
          return false;
      }
      if (dt->dt_stats.nr_differences > 5000 &&
          !nm.is_confluent() && !maf.options.force_differences)
      {
         /* if the number of word-differences is very high
            don't continue to favour equations with word-differences */
        stats.favour_differences = false;
        return false;
      }
    }
    return dt->is_interesting_equation(we,flags);
  }
  return false;
}

bool Rewriter_Machine::is_accepted_equation(const Working_Equation & we) const
{
  return dt->is_accepted_equation(we);
}

/**/

Rewriter_Machine::Rewriter_Machine(MAF & maf_) :
  pd(maf_.properties()),
  maf(maf_),
  container(maf_.container),
  nm(* new Node_Manager(*this)),
  em(* new Equation_Manager(*this)),
  jm(* new Job_Manager(*this)),
  nr_generators(maf_.generator_count()),
  nr_good_generators(maf_.generator_count()),
  nr_trivial_generators(0),
  nr_coset_reducible_generators(0),
  pool(0),
  derivation_db(0),
  dt(0),
  expand_ont(0),
  word_reducer(*new Equation_Word_Reducer(nm)),
  rm_state(RMS_Starting)
{
  memset(&stats,0,sizeof(stats));
  reset_pool();
  stats.first_timeout = maf.options.timeout/2;
  stats.normal_timeout = 2;
  stats.total_increment = 2;
  stats.priority_status = Priority_General;
  stats.want_overlaps = true;
  generator_properties = new unsigned char[nr_generators];
  memset(generator_properties,
         pd.has_cancellation ? GF_LEFT_CANCELS|GF_RIGHT_CANCELS : 0,
         nr_generators);
  /* We should set up the balancing flags before processing the axioms */
  if (maf.options.right_balancing == BYTE_OPTION_UNSET)
  {
    if (alphabet().order_properties[alphabet().order_type()].weight_moves_right())
      maf.options.right_balancing = 3;
    else
      maf.options.right_balancing = 1;
  }
  if (maf.options.left_balancing == BYTE_OPTION_UNSET)
  {
    if (maf.options.right_balancing == 0)
      maf.options.left_balancing = 0;
    else
      maf.options.left_balancing = alphabet().order_properties[alphabet().order_type()].weight_moves_right() ? 2 : 1;
  }
  maf.options.balancing_flags = 0;
  if (maf.options.right_balancing)
    maf.options.balancing_flags |= 1;
  if (maf.options.left_balancing)
    maf.options.balancing_flags |= 2;
  if (maf.options.no_equations)
    maf.options.no_deductions = true;
}

/**/

Rewriter_Machine::~Rewriter_Machine()
{
  Node_Reference e;
  while (expand_list.use(&e,nm))
    e->detach(nm,e);
  while (re_expand_list.use(&e,nm))
    e->detach(nm,e);

  if (dt)
  {
    delete dt;
    dt = 0;
  }
  if (expand_ont)
  {
    delete expand_ont;
    expand_ont = 0;
  }
  nm.destroy_tree();
  delete &jm;
  purge();
  delete &em;
  delete &nm;
  if (pool)
    delete pool;
  if (derivation_db)
    delete derivation_db;
  delete [] generator_properties;
  delete &word_reducer;
}

/**/

void Rewriter_Machine::start()
{
  stats.started = stats.starting = true;
  if (pd.is_group)
  {
    stats.want_differences = !maf.options.no_differences && pd.is_short ||
                             maf.options.force_differences || maf.options.differences;
  }
  else if (pd.is_coset_system)
  {
    if (maf.options.force_differences)
      stats.want_differences = true;
    else
      stats.want_differences = !maf.options.no_differences && pd.g_is_group &&
                               pd.g_level && pd.h_level || maf.options.differences;
  }
  else
  {
    maf.options.no_differences = true;
    stats.want_differences = false;
  }
  if (!pd.g_is_group)
  {
    maf.options.no_differences = true;
    maf.options.assume_confluent = true;
    stats.want_differences = false;
  }
  if (!pd.is_coset_system)
  {
    maf.options.no_h = false;
    if (maf.options.conjugation == 4)
      maf.options.conjugation = 1;
  }

  if (maf.options.assume_confluent)
    stats.want_differences = false;
  if (!stats.want_differences)
    maf.options.assume_confluent = true;
  if (pd.is_coset_system && maf.options.detect_finite_index)
    maf.options.assume_confluent = false;

  stats.want_secondaries = stats.want_differences;
  if (stats.want_secondaries)
    maf.options.secondary = alphabet().order_needs_moderation() ? 3 : 2;
  if (maf.options.secondary == 4)
    stats.want_secondaries = true;
  if (maf.options.secondary == BYTE_OPTION_UNSET)
    maf.options.secondary = alphabet().order_needs_moderation() ? 3 : 1;
  bool default_strategy = false;

  if (maf.options.strategy == ~0u)
  {
    default_strategy = true;
    maf.options.strategy = MSF_RESPECT_RIGHT;
    if (stats.want_differences && pd.g_level && !maf.options.tight)
      maf.options.strategy |= MSF_SELECTIVE_PROBE;
    if (!pd.inversion_difficult)
      maf.options.strategy |= MSF_AGGRESSIVE_DISCARD;
  }
  if (stats.want_differences)
  {
    if (maf.options.fast && maf.options.probe_style == BYTE_OPTION_UNSET)
      maf.options.probe_style = 0;
    if (!maf.options.tight && pd.g_level && maf.options.probe_style == BYTE_OPTION_UNSET)
      maf.options.probe_style = 2;
  }
  if (maf.options.probe_style == BYTE_OPTION_UNSET)
    maf.options.probe_style = 3;

  stats.favour_differences = stats.want_differences && !maf.options.tight;
  if (stats.want_differences)
  {
    if (!maf.options.is_kbprog && !alphabet().order_supports_automation())
      container.error_output("Warning! For this ordering MAF can compute word"
                             " differences, but automatic\nstructures"
                             " only if the rewriting system is confluent.\n");
    if (pd.inversion_difficult)
      maf.options.no_half_differences = true;
    if (!dt)
    {
      dt = Difference_Tracker::create(nm,1024);
      dt->close();
    }
  }
  else
    maf.options.no_kb = false;

  if (stats.visible_limit > 20 && !(maf.options.strategy & MSF_NO_FAVOUR_SHORT))
  {
    stats.visible_limit = stats.auto_expand_limit;
    /* If the spread of axiom sizes is not too great we will
       favour the shorter axioms by reducing the limits.
       However, this might not work well with long relators,
       if the axiom sizes are very different user should probably
       use era based expansion, or set -pool_above and -expand_all */
    if (stats.visible_limit > 2*stats.shortest_recent)
    {
      if (default_strategy)
        maf.options.strategy &= ~MSF_AGGRESSIVE_DISCARD;
    }

    if (stats.shortest_recent < stats.auto_expand_limit)
    {
      stats.auto_expand_limit = stats.shortest_recent;
      if (stats.visible_limit > stats.shortest_recent + 8)
        stats.visible_limit = stats.shortest_recent + 8;
    }
  }
//  if (maf.options.expand_all)
//    maf.options.strategy |= MSF_CLOSED;

  set_initial_special_limit();
  if (stats.visible_limit < maf.options.pool_above)
    stats.visible_limit = maf.options.pool_above;
  if (!stats.want_differences)
    maf.options.strategy |= MSF_CONJUGATE_FIRST;
  if (maf.options.strategy & MSF_CONJUGATE_FIRST)
    jm.divert_differences();
  if (maf.options.no_pool)
    stats.visible_limit = stats.auto_expand_limit = UNLIMITED;
  if (maf.options.no_deductions)
  {
    maf.options.no_half_differences = true;
    maf.options.check_inverses = false;
  }
  update_machine();
  stats.starting = false;
  if (maf.options.max_stored_length[0] &&
      stats.visible_limit < maf.options.max_stored_length[0]+maf.options.max_stored_length[1])
    stats.visible_limit = maf.options.max_stored_length[0]+maf.options.max_stored_length[1];
  stats.relaxed = !maf.options.tight || !stats.want_differences || !pd.g_level;
  stats.vital_attained = !stats.want_differences || !pd.g_level;
  if (stats.visible_limit > height()*2 && !maf.options.no_pool)
  {
    /* If a collapse has taken place already reduce the limits */
    stats.visible_limit = height()*2;
    if (stats.auto_expand_limit > stats.visible_limit)
      stats.auto_expand_limit = stats.visible_limit/2;
    if (stats.visible_limit < maf.options.pool_above)
      stats.visible_limit = maf.options.pool_above;
  }
  if (maf.options.strategy & MSF_NO_FAVOUR_SHORT)
    if (stats.visible_limit < stats.auto_expand_limit/5*8)
      stats.visible_limit = stats.auto_expand_limit/5*8;
  set_initial_special_limit();

  if (maf.options.expansion_order == 0)
    if (maf.options.strategy & MSF_USE_ERAS)
      maf.options.expansion_order = 4;
    else if (pd.is_short || !pd.inversion_difficult)
      maf.options.expansion_order = 1;
    else
      maf.options.expansion_order = 2;

  switch (maf.options.expansion_order)
  {
    case 1:
      break;
    case 2:
      expand_ont = new By_Size_Node_Tree(nm);
      break;
    case 3:
      expand_ont = new By_SizeTime_Node_Tree(nm);
      break;
    case 4:
      expand_ont = new By_Time_Node_Tree(nm);
      break;
    case 5:
    case 6:
      expand_ont = new By_Left_Size_Node_Tree(nm,maf.options.expansion_order==6);
      break;
  }
  rm_state = RMS_Starting;
}

void Rewriter_Machine::reset_pool()
{
  pool = 0;
  stats.pool_limit = stats.shortest_recent = UNLIMITED;
  stats.best_pool_lhs = stats.best_pool_rhs = WHOLE_WORD;
  stats.best_pool_lhs_total = stats.best_pool_rhs_total = UNLIMITED;
  stats.discard_limit = stats.pool_limit + POOL_DEPTH;
}

/**/

void Rewriter_Machine::restart()
{
  /* This method is called if update_machine() believes a "collapse" has
     taken place, and that the contents of the pool are worthless, or if
     the pool has become too big.
     In these cases we throw away the pool completely, and restart KB from
     scratch. We have to add the axioms again in case we would otherwise
     have thrown away something necessary.
     This method can save a lot of time in the confluent case, especially
     if a recursive ordering is used. */
  if (pool)
  {
    delete pool;
    reset_pool();
  }
  Node_Reference e;
  stats.visible_limit = 0;
  stats.auto_expand_limit = 0;
  stats.lhs_limit = 0;
  Node_Iterator ni(nm);
  ni.scan(&e);
  while (ni.next(&e))
    if (e->is_final() && e->fast_is_primary())
    {
      if (e->status() == NS_Expanded_Equation)
      {
        Total_Length tl = e->total_length(nm);
        if (tl > stats.visible_limit)
          stats.visible_limit = tl;
        if (tl > stats.auto_expand_limit)
          stats.auto_expand_limit = tl;
        Word_Length l = e->length(nm);
        if (l > stats.lhs_limit)
          stats.lhs_limit = l;
        e->node(nm).reduction.expand_timestamp = 0;
      }
    }
  stats.max_expanded_lhs = stats.lhs_limit;
  stats.special_limit = stats.max_expanded = stats.auto_expand_limit;
  stats.approved_complete = stats.complete = stats.primary_complete = false;
  const Linked_Packed_Equation *axiom;
  for (axiom = maf.first_axiom();axiom;axiom = axiom->get_next())
  {
    Simple_Equation se(alphabet(),*axiom);
    add_axiom(se.lhs_word,se.rhs_word,0);
  }
  set_initial_special_limit();
}

/**/

bool Rewriter_Machine::short_status(int priority,unsigned level,int gap,const char * control,Variadic_Arguments & va)
{
  if (stats.priority_status > priority)
    gap = 0;

  if (container.status(level,gap,control,va))
  {
    stats.status_count += gap;
    stats.priority_status = priority;
    return true;
  }
  else if (container.status_needed(gap))
  {
    stats.status_count += gap;
    return true;
  }
  return false;
}

/**/

bool Rewriter_Machine::container_status(unsigned level,int gap,const char * control,...)
{
  DECLARE_VA(va,control);
  return short_status(Priority_General,level,gap,control,va);
}

/**/

bool Rewriter_Machine::critical_status(unsigned level,int gap,const char * control,...)
{
  DECLARE_VA(va,control);
  return short_status(Priority_Critical,level,gap,control,va);
}

/**/

bool Rewriter_Machine::status(int priority,unsigned level,int gap,const char * control,...)
{
  DECLARE_VA(va,control);
  if (stats.priority_status > priority/* && stats.priority_status != 0*/)
    gap = 0;

  if (container.status(level,gap,control,va))
  {
    const Equation_Stats & estats = jm.status();
    stats.priority_status = priority;
    container.progress(2,"Nodes:L0=" FMT_NC ",L1=" FMT_NC ",L2/3=" FMT_NC ",Bad=" FMT_NC ",Depth=%d,Acc=%d\n",
                       nm.stats.nc[language_L0],nm.stats.nc[language_L1],
                       nm.stats.nc[language_L2]+nm.stats.nc[language_L3],
                       nm.stats.nc[language_A],height(),
                       nm.start()->reduced.max_accepted);
    if (stats.want_differences && dt)
    {
      const Difference_Tracker::Status & dstats = dt->status();
      long t = stats.timeout-(stats.status_count-stats.ticked);
      if (t < 0)
        t = 0;
      container.progress(2,"Timeout:%ld/%lu/%lu/%lu "
                           "Deltas:E=" FMT_NC ",D=" FMT_ID " Diffs:" FMT_ID "/" FMT_ID
                           " Sizes:%d/%d/%d/%d\n",
                         t,
                         stats.timeout,
                         stats.status_count - stats.last_difference,
                         stats.status_count - stats.update_time,
                         stats.equation_delta,stats.difference_delta,
                         dstats.nr_primary_differences,
                         dstats.nr_differences,
                         dstats.max_difference_length,
                         dstats.max_difference_eqn,
                         dstats.max_primary_difference_eqn,
                         dstats.max_interesting_eqn);
    }
    if (pd.is_coset_system)
      if (pd.nr_generators > pd.coset_symbol+1)
        container.progress(2,"Comp: G=%d,_H=%d,H=%d ",stats.g_complete,stats.coset_complete,stats.h_complete);
      else
        container.progress(2,"Comp: G=%d,_H=%d ",stats.g_complete,stats.coset_complete);

    if (!maf.options.no_pool)
      container.progress(2,"Limits:E=%d/%d,S/T/P/D=%d/%d/%d/%d,L=%d/%d/%d,O=%d,C=%d,F=%d\n",
                           stats.auto_expand_limit,stats.max_expanded,
                           stats.special_limit,stats.visible_limit,
                           stats.pool_limit,stats.discard_limit,
                           stats.lhs_limit,
                           stats.max_expanded_lhs,
                           stats.max_primary_lhs,
                           nm.overlap_filter.overlap_limit,
                           stats.upper_consider_limit,
                           nm.overlap_filter.forbidden_limit);
    if (pool)
      container.progress(2,"Pool: Best LHS %d/%d Best RHS %d/%d\n",
                         stats.best_pool_lhs,stats.best_pool_lhs_total,
                         stats.best_pool_rhs,stats.best_pool_rhs_total);
    container.progress(2,"Equations:Total=" FMT_NC ", Using=" FMT_NC
                         ", Visible=" FMT_NC ", Pool=" FMT_ID ".\n",
                       estats.nr_equations,estats.nr_adopted,estats.nr_visible,
                       pool ? pool->count() : 0);

    stats.status_count += gap;
    return true;
  }
  else if (container.status_needed(gap))
  {
    stats.status_count += gap;
    return true;
  }
  return false;
}

void Rewriter_Machine::late_differences()
{
  /* If we were in assume_confluent mode, then we did not calculate word
     differences, and did not create all the conjugate equations.
     So we do both of those now */

  if (!maf.options.no_differences)
    stats.favour_differences = stats.want_differences = true;
  if (pd.inversion_difficult && !stats.complete)
    maf.options.no_half_differences = true;
  if (stats.want_differences)
  {
    stats.want_secondaries = true;
    if (stats.complete)
      stats.want_overlaps = false;
    maf.options.secondary = pd.inversion_difficult ? 3 : 2;
    if (!dt)
    {
      dt = Difference_Tracker::create(nm,1024);
      dt->close(true);


      Node_Reference e;
      Node_Count to_do = jm.status().nr_visible;
      Node_Count done = 1;
      int gap = 2;
      Node_Iterator ni(nm);
      ni.scan(&e);
      while (ni.next(&e))
      {
        if (e->is_final() && (!pd.is_coset_system || e->last_letter() < pd.coset_symbol))
        {
          done++;
          schedule_right_conjugation(e);
          coset_extras(e);
          if (!(char) done && container_status(2,gap,
                      "Conjugating equations (" FMT_NC " of " FMT_NC ")\n",++done,to_do))
            gap = 1;
        }
      }
      update_machine();
    }
    check_differences(CD_UPDATE);
  }
}

/**/

int Rewriter_Machine::expand_machine(bool limit)
{
  stats.visibility_correct = false; // clear this flag, which might be set if
                                    // we are resuming expansion after a failed
                                    // axiom check
  /* the next three lines are needed in case we are resuming expansion of a coset system
     for a finite index subgroup and -detect_finite_index was not originally specified */
  stats.no_coset_pool = pd.is_coset_finite;
  if (stats.no_coset_pool)
    rm_state |= RMS_Stick;

  update_machine(limit ? UM_CHECK_POOL : 0);
  if (limit)
  {
    stats.approved_complete = false;
    maf.options.strategy |= MSF_AGGRESSIVE_DISCARD;
  }
  if (maf.options.no_kb)
  {
    maf.options.no_kb = false;
    Total_Length vital_limit = dt->vital_eqn_limit();
    if (stats.visible_limit < vital_limit)
      stats.visible_limit = vital_limit;
    set_initial_special_limit();
    return 1;
  }
  if (!limit)
    container.progress(1,"Looking for new equations using Knuth Bendix process\n");
  for (;;)
  {
    Node_Reference e;

    schedule(limit ? rm_state|RMS_Limited : rm_state);
    Node_Count done = 0;
    stats.ticked = stats.status_count;
    stats.eticked = stats.status_count;
    stats.expand_time = stats.status_count+1;
    stats.timeout = stats.difference_delta &&
              stats.status_count-stats.last_difference < maf.options.timeout ?
                    stats.normal_timeout : stats.first_timeout;
    rm_state = 0;
    stats.shortest_recent = UNLIMITED;
    stats.just_ticked = false;
    stats.hidden_differences = false;

    int phase = expand_ont ? 0 : 1;
    stats.start_id = nm.last_id;
    stats.some_discarded = false;
    Overlap_Search_Type ost = /*limit ? OST_Quick :*/ OST_Full;

    for (;;)
    {
      bool found = false;
      if (expand_ont)
        found = expand_ont->use(&e);
      else
      {
        if (phase == 1)
        {
          found = expand_list.use(&e,nm);
          if (!found)
          {
            if (!(rm_state & RMS_Abnormal))
            {
              if (re_expand_list.length() > 1000 && pool &&
                  height() > stats.visible_limit/2 + 2 &&
                  dt && dt->status().nr_differences < 5000)
              {
                update_machine(UM_CHECK_POOL+UM_PRUNE);
                dt->scan_all_differences();
              }
            }
            if (rm_state == RMS_Hurrying)
            {
              rm_state = RMS_Normal;
              ost = OST_Full;
            }
            phase = 2;
          }
        }
        if (phase == 2)
          found = re_expand_list.use(&e,nm);
      }
      if (!found)
        break;
      ++done;
      e->node(nm).clear_flags(EQ_EXPANDED);
      bool is_coset_equation = false;
      if (e->status() == NS_Expanded_Equation)
      {
        if (maf.aborting)
          rm_state = RMS_Aborting;

        if (!(rm_state & RMS_Aborting))
        {
          is_coset_equation = pd.is_coset_system &&
                              e->first_letter(nm) >= pd.coset_symbol &&
                              e->last_letter() < pd.coset_symbol;
          expand(e,done,
                 is_coset_equation && !e->expanded_timestamp() ? OST_Full : ost);
          if (maf.options.repeat == 1 && e->status() == NS_Expanded_Equation)
            expand(e,done,ost);
        }

        if (stats.status_count-stats.eticked > 2)
        {
          stats.shortest_recent = UNLIMITED;
          if (!stats.complete && !stats.equation_delta &&
              !stats.no_coset_pool && !limit && rm_state == RMS_Normal &&
              (nm.has_filtered() || stats.some_discarded))
          {
            ost = OST_Aborted_Pass;
            rm_state |= RMS_Aborted_Pass|RMS_Widen;
            container.progress(2,"Increasing limits since no new equations recently\n");
          }
          stats.eticked = stats.status_count;
          stats.equation_delta /= 2;
          if (pd.is_coset_system && maf.options.detect_finite_index)
          {
            if (stats.status_count - stats.last_coset_equation > 10)
            {
              stats.coset_complete = true;
              if (rm_state == RMS_Normal)
              {
                if (is_coset_equation && !stats.no_coset_pool &&
                    stats.status_count - stats.last_coset_equation > 20 ||
                    nr_coset_reducible_generators == pd.coset_symbol)
                {
                  ost = OST_Aborted_Pass;
                  rm_state |= RMS_Aborted_Pass|RMS_Want_Finite_Check;
                }
              }
            }
            else if (!stats.coset_complete && rm_state & RMS_Want_Finite_Check)
            {
              ost = OST_Full;
              rm_state = RMS_Normal;
            }
          }
          if (dt && (
              maf.options.max_equations && nm.stats.nc[1] >= maf.options.max_equations ||
              maf.options.max_time && stats.status_count >= maf.options.max_time))
          {
            ost = OST_Aborted_Pass;
            rm_state |= RMS_Aborted_Pass;
          }
        }
        if (dt && !maf.options.assume_confluent && !(rm_state & RMS_Abnormal))
        {
          if (!maf.options.fast && pd.g_level && stats.favour_differences &&
              dt->vital_eqn_limit() > stats.visible_limit &&
              e->total_length(nm) > dt->vital_eqn_limit())
          {
            /* The maximum word-difference length has increased past the size
               of equations we are keeping. So we need to increase
               stats.visible_limit to ensure we get all the transitions.
            */
            stats.vital_attained = false;
            rm_state |= RMS_Too_Tight;
            ost = OST_Aborted_Pass;
          }
          if (!is_dm_changing())
          {
            if (stats.just_ticked)
            {
              stats.just_ticked = false;
              switch (ost)
              {
                case OST_Full:
                  if (stats.complete || !pool || maf.options.slow)
                    break;
                  ost = OST_Moderate;
                  break;
                case OST_Moderate:
                  ost = OST_Quick;
                  break;
                case OST_Quick:
                  if ((phase==2 ||
                       stats.status_count-stats.last_difference >= maf.options.timeout &&
                       stats.status_count-stats.expand_time >= maf.options.timeout) &&
                      build_allowed(false))
                  {
                    ost = OST_Minimal;
                    if (phase == 2)
                      container.progress(1,"Will attempt to build automata at end of pass\n");
                  }
                  break;
                case OST_Minimal:
                case OST_Aborted_Pass:
                  break;
              }
              if (ost != OST_Full)
                rm_state |= RMS_Hurrying;
            }
          }
          else
          {
            switch (ost)
            {
              case OST_Moderate:
                if (stats.just_ticked)
                  ost = OST_Full;
                break;
              case OST_Quick:
                if (stats.just_ticked)
                  ost = OST_Moderate;
                break;
              case OST_Minimal:
                ost = OST_Quick;
                break;
              case OST_Full:
              case OST_Aborted_Pass:
                break;
            }
            stats.just_ticked = false;
            if (phase==1 || ost == OST_Full)
              rm_state &= ~RMS_Hurrying;
          }
        }
      }
      e->detach(nm,e);
    }
    if (nm.has_filtered())
      stats.complete = false;
    bool should_build = examine(rm_state);
    bool coset_complete_flag = should_build && stats.no_coset_pool;

    if (!should_build && pd.is_coset_system && stats.coset_complete &&
        !stats.complete &&
        maf.options.detect_finite_index && !maf.options.assume_confluent)
    {
      container_status(2,1,"Checking current size of coset language\n");
      FSA_Simple * temp = 0;
      Rewriting_System * temp_rws = new
         Rewriting_System(this,RWSC_MINIMAL|RWSC_FSA_ONLY|RWSC_NEED_FSA|RWSC_NO_H,&temp);
      delete temp_rws;
      FSA_Simple * group_wa = FSA_Factory::restriction(*temp,maf.group_alphabet());
      Word_Length max_coset_length;
      Language_Size count = group_wa->language_size(&max_coset_length,false);
      delete group_wa;
      if (count <= 2)
      {
        delete temp;
        maf.options.assume_confluent = true;
        maf.options.no_h = false;
        if (stats.auto_expand_limit <= max_coset_length*2+4)
          stats.auto_expand_limit = max_coset_length*2+4;
        should_build = stats.no_coset_pool = false;
        container.progress(1,"Group is finite - will continue until confluence\n");
      }
      else
      {
        temp->set_is_initial(1,false); // do this first in case the subgroup has no generators
        temp->set_is_initial(temp->new_state(1,pd.coset_symbol),true);
        // turn group_wa into a coset word-acceptor
        group_wa = FSA_Factory::restriction(*temp,maf.group_alphabet());
        delete temp;
        count = group_wa->language_size(&max_coset_length,true);
        if (count != LS_INFINITE)
        {
          should_build = stats.no_coset_pool;
          stats.no_coset_pool = coset_complete_flag = true;
          if (stats.auto_expand_limit <= max_coset_length*2+4)
            stats.auto_expand_limit = max_coset_length*2+4;
          container.progress(1,"Index is finite (" FMT_LS ")\n",count);
          hide_equations(group_wa);
          if (should_build)
          {
            update_machine(UM_CHECK_POOL);
            late_differences();
            stats.coset_complete = true;
          }
        }
        else
          hide_equations(group_wa);
      }
    }
    if (stats.complete)
    {
      nm.validate_tree();
      if (!limit)
        container.progress(1,"All possible primary equations have been deduced\n");
      nm.confluent = true;
      late_differences();
      stats.favour_differences = stats.want_differences;
      purge();
      return 0;
    }

    if (should_build || limit && stats.approved_complete)
    {
      stats.coset_complete = coset_complete_flag;
      stats.favour_differences = stats.want_differences;
      return maf.aborting ? -1 : 1; /* We will try to build the automatic structure */
    }
    else if (stats.status_count - stats.last_save_time > 120 &&
             !stats.no_coset_pool && !maf.options.emulate_kbmag)
    {
      stats.last_save_time = stats.status_count;
      return 2;
    }
    /* we should do this here rather than in the test for finite coset system
       because if a new coset equation is discovered after this is set then
       stats.coset_complete might not be set any more, but we still do not
       want to increment the limits */
    if (pd.is_coset_system && stats.no_coset_pool)
      rm_state |= RMS_Stick;
  }
}

/**/

void Rewriter_Machine::schedule(unsigned rm_state)
{
  /* This method adjusts the various limits that are currently in effect,
     and then builds a new list of equations to be KB expanded, possibly
     after examining the pool first.
     This method is rather messy as there are a lot of limits to be adjusted
     and we have to be very careful what we do, otherwise the program can
     start misbehaving */

  Total_Length old_visible_limit = stats.visible_limit;
  bool increment_expand_limits = stats.approved_complete ||
                                 rm_state & RMS_Widen ||
                                 stats.relaxed && !(rm_state & RMS_Starting+RMS_Too_Tight);
  if (rm_state & RMS_Stick+RMS_Limited)
    increment_expand_limits = false;
  Total_Length vital_eqn_limit = 0;
  Total_Length upper_consider_limit = 0;
  stats.pool_full = false;
  bool pool_checked = false;
  Element_Count pool_count = pool ? pool->count() : 0;
  if (!(rm_state & RMS_Limited+RMS_Stick))
  {
    if (dt && !maf.options.assume_confluent && pd.g_level && stats.favour_differences)
    {
      /* In this case we are presumably looking for a shortlex type automatic
         structure. We can't use this code to set the limits if we are looking
         for some other kind of automatic structure, because the length
         of a word-difference is unrelated to the size of the equation it
         comes from, especially if we are using a recursive type ordering

         Here we set the size up to which equations are automatically inserted
         in the Rewriter_Machine. Frequently most word-differences arise in
         fairly short equations, but a few are only in very long equations.
         Until we can build the automata successfully we can't be sure that
         apparent word-differences aren't spurious ones caused by our missing
         some shorter unknown word-difference state or transition.

         If we increase this limit either too slowly or too fast we will tend
         to get more spurious differences, in the former case because we lack
         enough equations even to reduce only moderately long words correctly,
         in the latter case because we are constructing equations in words that
         are very long. Also, in the former case our tree of equations will be
         unnecessarily thick, and in the latter case it will have excessively
         long branches that will impact the rate at which new equations can be
         inserted. We want to try to avoid letting the visible_limit get
         too big if the object is infinite, and to avoid stopping it getting
         big enough if the object is finite!
      */

      vital_eqn_limit = dt->vital_eqn_limit();
      /* Make sure we expand all equations that might impact the computation
         of the transitions in the word-difference machine */
      if (stats.auto_expand_limit < vital_eqn_limit && !increment_expand_limits)
      {
        is_dm_changing();
        check_differences(CD_CHECK_INTEREST);
        vital_eqn_limit = dt->vital_eqn_limit();
      }

      if (stats.auto_expand_limit < vital_eqn_limit)
      {
        if (stats.auto_expand_limit * 2 < stats.visible_limit)
        {
          stats.auto_expand_limit *= 2;
          increment_expand_limits = false;
        }
        else if (stats.auto_expand_limit + 4 < vital_eqn_limit)
        {
          stats.auto_expand_limit += 4;
          increment_expand_limits = false;
        }
        else if (stats.visible_limit >= vital_eqn_limit)
        {
          stats.auto_expand_limit = vital_eqn_limit;
          increment_expand_limits = false;
        }
        else
        {
          /* In this case we cannot increase the auto_expand_limit safely.
             We are likely to get here if the group is not really automatic */
          if (stats.auto_expand_limit < stats.visible_limit)
          {
            stats.auto_expand_limit = stats.visible_limit;
            increment_expand_limits = false;
          }
        }
      }

      /* We increase visible_limit towards the interest limit, at a rate that
         depends on how far below the limits we are */
      const Difference_Tracker::Status &dstats = dt->status();
      if (stats.visible_limit < stats.auto_expand_limit*2)
      {
        if (!(rm_state & RMS_Too_Tight))
        {
          if (stats.visible_limit < dstats.max_primary_difference_eqn)
            stats.visible_limit++;
          if (stats.visible_limit < dstats.max_difference_eqn)
            stats.visible_limit++;
          if (stats.visible_limit < dstats.max_interesting_eqn)
            stats.visible_limit++;
        }
      }
      if (upper_consider_limit < dstats.max_primary_difference_eqn)
        upper_consider_limit = dstats.max_primary_difference_eqn;
      if (upper_consider_limit < dstats.max_difference_eqn)
        upper_consider_limit = dstats.max_primary_difference_eqn;

      if (rm_state & RMS_Too_Tight)
      {
        /* We won't let height limit stay below the vital size */
        if (stats.visible_limit < vital_eqn_limit)
          stats.visible_limit = vital_eqn_limit;
      }

      /* don't let visible limit have a silly value */
      if (stats.visible_limit & 1 && stats.total_increment==2)
        stats.visible_limit++;
      /* The tests below are for the situation where it looks like we
         are going to get confluence. In this case we want to try to
         stop KB being aborted */
      /* If there is no pool we may as well expand all the equations */
      if (!(rm_state & RMS_Too_Tight) && increment_expand_limits)
        if (!pool && stats.auto_expand_limit < stats.max_expanded)
          stats.auto_expand_limit = stats.max_expanded;
      /* If the pool is very small, and we have a fair number of
         equations already try to get all of it altogether, since
         this will stop early attempts to build */
      Total_Length upper_reasonable_limit = min(stats.discard_limit,stats.upper_consider_limit);
      if (pool && pool->count() < 100 && nm.stats.nc[1] > 200 &&
          stats.visible_limit < upper_reasonable_limit &&
          !(rm_state & RMS_Too_Tight+RMS_Starting))
        stats.visible_limit = upper_reasonable_limit;
    }
    else if (nr_good_generators > 1 && pool)
    {
      Element_Count count = pool_count;
      while (count/nr_good_generators > 100)
      {
        count /= nr_good_generators;
        stats.visible_limit += stats.total_increment;
      }

      if (increment_expand_limits && pd.g_level)
      {
        /* In shortlex case, and especially while tree is small, increase the
           limits a bit more quickly */
        stats.visible_limit += stats.total_increment;
        if (nm.stats.nc[language_L1] < 10000)
          stats.visible_limit += stats.total_increment;
        if (nm.stats.nc[language_L1] < 1000)
          stats.visible_limit += stats.total_increment;
      }
    }
    else
      stats.visible_limit += stats.total_increment;
  }

  /* detect the unusual situation where few if any equations were being
     constructed due to limits being overly restrictive. We increase the
     height limit to take on more of the pool */
  if (rm_state & RMS_Widen && !stats.some_discarded &&
      stats.visible_limit < stats.upper_consider_limit)
    stats.visible_limit = stats.upper_consider_limit;

  if (increment_expand_limits)
  {
    /* Push the lower limits towards the upper ones if need be */
    if (stats.approved_complete && (!stats.want_differences || !pd.g_level))
    {
      if (stats.auto_expand_limit < stats.visible_limit)
        stats.auto_expand_limit += stats.total_increment;
      if (stats.auto_expand_limit < stats.max_expanded)
        stats.auto_expand_limit += stats.total_increment;
    }

    stats.auto_expand_limit += stats.total_increment;

    if (stats.max_expanded < stats.auto_expand_limit)
      stats.max_expanded = stats.auto_expand_limit; /* to avoid confusing pool check */
    if (rm_state & RMS_Widen)
    {
      /* relax overly restrictive limits */
      Total_Length new_limit = min(stats.max_expanded,stats.visible_limit);
      if (stats.auto_expand_limit < new_limit)
        stats.auto_expand_limit = new_limit;
    }
  }
  Word_Length shortest_fresh_lhs = WHOLE_WORD;
  for (;;)
  {
    /* Keep the height limit away from the expand_limit */
    if (!(rm_state & RMS_Limited))
    {
      if (stats.visible_limit < stats.auto_expand_limit+4)
        stats.visible_limit = stats.auto_expand_limit+4;
      /* detect the situation where a generator should be eliminated but the
         rhs for it is longish */
      if (stats.best_pool_lhs == 1 && stats.best_pool_lhs_total < stats.discard_limit)
        if (stats.visible_limit < stats.best_pool_lhs_total)
          stats.visible_limit = stats.best_pool_lhs_total;
    }

    if (stats.pool_limit <= stats.visible_limit)
    {
      /* We think there is something in the pool worth having */
      old_visible_limit = stats.visible_limit;
      set_initial_special_limit();
      update_machine(UM_CHECK_POOL+UM_PRUNE);
      pool_checked = true;
      if (pool && dt &&
          stats.favour_differences &&
          nm.stats.nc[0] < 500000 && dt->dm_changed(0) &&
          dt->status().nr_differences < 4000 &&
          !(rm_state & RMS_Starting+RMS_Too_Tight) &&
          stats.last_explore + stats.explore_time*2 < stats.status_count)
      {
        /* In this case our main tree is unreliable, because we have
           very likely not being doing the KB expansion methodically.
           So we use explore_dm() to kill off reducible nodes, and hopefully
           bad differences. Then, since we are favouring equations with
           differences make sure that these have been calculated appropriately*/
        explore_dm(false);
        check_differences(CD_CHECK_INTEREST);
      }
      if (dt && vital_eqn_limit)
        vital_eqn_limit = dt->vital_eqn_limit();
    }
    stats.complete = !pool;
    /* some_ignored is used only to track discards that were forced by
       user limits, not our own limits */
    stats.some_ignored = false;

    /* Now we have some provisional limits, and take a look at what
       equations we have. We build a list of all the primaries and
       choose some to expand */
    stats.approved_complete = true;
    stats.g_complete = stats.h_complete = stats.coset_complete = true;
    Node_Reference e;
    Multi_Node_List mnl_new(nm,height()),mnl_old(nm,height());
    Node_List nl;
    Node_List onl;
    Node_Count to_do = primary_count();
    Node_Count done = 0;
    Total_Length non_expanded_limit = UNLIMITED;
    int gap = 2;
    Node_Iterator ni(nm);
    while (ni.scan(&e))
    {
      if (e->is_equation())
      {
        if (e->status() == NS_Oversized_Equation)
        {
          jm.cancel_jobs(e);
          jm.schedule_jobs(e,choose_status(e,NS_Undifferenced_Equation));
        }
        Total_Length tl = e->total_length(nm);
        if (tl <= stats.visible_limit && tl > old_visible_limit)
          schedule_right_conjugation(e);
        if (e->fast_is_primary())
        {
          if (container_status(2,gap,"Sorting primary equations (" FMT_NC
                                     " of " FMT_NC ")\n",++done,to_do))
            gap = 1;
          if (!e->expanded_timestamp() && tl < non_expanded_limit)
            non_expanded_limit = tl;
          mnl_new.add_state(e);
        }
        else
        {
          /* We schedule a check of all the secondaries again */
          if (e->raw_reduced_node(nm)->is_final())
            nm.publish_bad_rhs(e);
          else
            nm.publish_dubious_rhs(e);
        }
      }
    }
    mnl_new.collapse(&nl);
    if (!maf.options.assume_confluent)
      if (non_expanded_limit > stats.auto_expand_limit+stats.total_increment &&
          non_expanded_limit < stats.max_expanded)
      {
        /* In this case we find that we have in fact expanded all the
           equations up to a size above the expand_limit. So we increase it.
           This may help the limits grow quicker when we are close to
           confluence.*/
        stats.auto_expand_limit = non_expanded_limit - stats.total_increment;
      }
    bool force = maf.options.expand_all;
    Node_Count nr_fresh = 0;
    Node_Count nr_scheduled = 0;
    stats.max_expanded_lhs = 0;
    stats.lhs_limit = 0;
    stats.max_primary_lhs = 0;
    // It seems to be bad to reset max_expanded.
    //stats.max_expanded = 0;

    bool again = true;
    for (;;)
    {
      to_do = nl.length();
      done = 0;
      Total_Length non_expand_limit = UNLIMITED;
      stats.vital_attained = stats.visible_limit >= vital_eqn_limit;
      while (nl.use(&e,nm))
      {
        if (!(char ) done)
          container_status(2,1,"Scheduling equations ("
                               FMT_NC " of " FMT_NC ")\n",++done,to_do);
        Word_Length l = e->lhs_length(nm);
        Total_Length tl = e->total_length(nm);
        if (l > stats.max_primary_lhs)
          stats.max_primary_lhs = l;
        State rhs = e->raw_reduced_node(nm);

        if (rhs->is_final())
          nm.publish_bad_rhs(e); /* Have another bash at repairing a broken equation */
        else if (e->status() <= NS_Adopted_Equation &&
                 schedule_expand(mnl_new,mnl_old,e,force))
        {
          if (e->expanded_timestamp()==0)
          {
            nr_fresh++;
            if (l < shortest_fresh_lhs)
              shortest_fresh_lhs = l;
          }
          nr_scheduled++;
        }
        else
        {
          if (tl < non_expand_limit)
            non_expand_limit = tl;
          onl.add_state(e);
        }
      }
      nl.merge(&onl);
      to_do = nl.length();
      if (to_do)
      {
        /* Here we find that some equations are not going to be expanded */
        if (again)
        {
          /* First time round lhs_limit was not properly known, so we may
             have missed some equations we would like to expand in the
             non-shortlex case. So we go round again to pick these up */
          again = false;
          continue;
        }
        /* If there are very few equations left that have never been expanded,
           just expand them all */
        if ((to_do < nr_scheduled/4 || to_do < 100) && !pool &&
            !(rm_state & RMS_Stick))
        {
          force = true;
          continue;
        }

        if (nr_fresh==0 && increment_expand_limits)
          if (non_expand_limit < stats.pool_limit)
          {
            /* Here we find that we have no new equations to expand. So
               we force the expand limit up and go around again, but only
               if the non-expanded equations are smaller than what is in the
               pool
            */
            stats.auto_expand_limit = non_expand_limit;
            again = true;
            continue;
          }
          else
            stats.auto_expand_limit = stats.pool_limit;
        stats.complete = false;
      }
      break;
    }
    if (!expand_ont)
    {
      mnl_new.collapse(&expand_list);
      mnl_old.collapse(&re_expand_list);
    }
    if (nr_fresh || !pool && stats.pool_limit == UNLIMITED || !increment_expand_limits)
    {
      /* We are happy with our list.

         Check the list of equations we aren't going to expand this time.
         If any were previously expanded we need to clear their timestamp */
      while (nl.use(&e,nm))
        if (e->status() == NS_Expanded_Equation)
          e->node(nm).reduction.expand_timestamp = 0;
      break;
    }
    /* We are not happy with the list - there are no new equations on it. So
       we throw it away, and are going to look in the pool for new equations */
    while (re_expand_list.use(&e,nm))
    {
      e->node(nm).clear_flags(EQ_EXPANDED);
      e->detach(nm,e);
    }
    if (expand_ont)
      while (expand_ont->use(&e))
      {
        e->node(nm).clear_flags(EQ_EXPANDED);
        e->detach(nm,e);
      }
    if (stats.visible_limit < stats.pool_limit)
      stats.visible_limit = stats.pool_limit;
  }

  /* We now have a list of equations to expand, which are possibly
     of very different sizes. We want to limit things so that
     we do not start creating lots of very long equations before we
     have found all the short ones.

     The values passed as parameters to nm.set_overlap_limits()
     can make a huge difference to MAF in the case where we are looking
     for a confluent RWS. The neql value especially influential
     even though it is decisive in only a small number of overlaps.

     We will at least consider all overlaps between LHSes up to an
     overlap of 1 between what we would expect for ths LHS length that
     occurs in equations <= the auto_expand_limit. We will most likely
     be expanding some longer equations, and we ensure that we go at least
     two beyond the longest such LHS.

     We also expand overlaps that give an equation (prior to reduction)
     up to a reasonable size. In the non-shortlex case some of these may
     come from unusually long overlaps.

     In the case where we are not assuming confluence, we may choose to
     expand longer overlaps. But we set a limit on this - overlaps cannot
     be longer than the longest possible overlap that existed at the start
     of the phase. Without such a limit, we would tend to create lots of
     long equations with possibly spurious word-differences. Also each phase
     would take a lot longer.
  */

  Total_Length lower_overlap_limit = stats.auto_expand_limit;
  Total_Length consider_limit;
  Total_Length upper_overlap_limit = stats.max_expanded_lhs+stats.max_primary_lhs;

  if (lower_overlap_limit < stats.max_expanded_lhs + 2)
    lower_overlap_limit = stats.max_expanded_lhs + 2;
  if (shortest_fresh_lhs != WHOLE_WORD &&
      lower_overlap_limit < stats.lhs_limit + shortest_fresh_lhs-1)
    lower_overlap_limit = stats.lhs_limit + shortest_fresh_lhs-1;
  /* We now finally adjust the height limits to keep them above
     what would get expanded */
  Total_Length nl = min(stats.lhs_limit + stats.auto_expand_limit/2,
                        stats.max_expanded+stats.total_increment);
  if (stats.visible_limit < nl)
    stats.visible_limit = nl;
  if (stats.discard_limit < stats.visible_limit + POOL_DEPTH)
    stats.discard_limit = stats.visible_limit + POOL_DEPTH;
  consider_limit = stats.visible_limit + POOL_DEPTH;
  if (stats.want_differences && pd.g_level &&
      stats.max_expanded+2 > consider_limit)
    consider_limit = stats.max_expanded+2;
  if (consider_limit < upper_consider_limit)
    consider_limit = upper_consider_limit;
  stats.upper_consider_limit = consider_limit+2;
  if (!(maf.options.filters & 1))
    lower_overlap_limit = upper_overlap_limit = UNLIMITED;
  if (!(maf.options.filters & 2))
    stats.upper_consider_limit = UNLIMITED;
  nm.set_overlap_limits(lower_overlap_limit,
                        stats.upper_consider_limit,
                        stats.visible_limit, upper_overlap_limit,
                        maf.options.probe_style & 1 ?
                        stats.max_primary_lhs : stats.max_expanded_lhs);
  for (Ordinal g = 0; g < nr_generators;g++)
    generator_properties[g] &= ~GF_DIFFICULT;
  set_initial_special_limit();
  update_machine(stats.pool_limit <= stats.visible_limit || !pool_checked ? UM_CHECK_POOL : 0);
}

/**/

void Rewriter_Machine::hide_equations(FSA_Simple * coset_word_acceptor)
{
  /* In a coset system all the g word equations for the group will be
     created by KB expansion. But some, perhaps most, will never be used
     for coset reduction. We look for such equations by finding nodes
     which correspond to words that cannot occur as part of a coset word,
     and set a flag which the code that builds word-differences checks
     and uses to ignore irrelevant equations.
  */

  if (pd.is_coset_system && !stats.visibility_correct)
  {
    if (!coset_word_acceptor)
    {
      FSA_Simple * temp = 0;
      Rewriting_System * temp_rws = new Rewriting_System(this,RWSC_FSA_ONLY|RWSC_NEED_FSA|RWSC_MINIMAL|RWSC_NO_H,&temp);
      delete temp_rws;
      temp->set_is_initial(1,false); // do this first in case the subgroup has no generators
      temp->set_is_initial(temp->new_state(1,pd.coset_symbol),true);
      coset_word_acceptor = FSA_Factory::restriction(*temp,maf.group_alphabet());
      delete temp;
    }
    if (coset_word_acceptor->language_size(false)==2 || nm.is_confluent())
      stats.visibility_correct = true;
    coset_word_acceptor->set_initial_all();
    FSA_Simple * visible_words = FSA_Factory::determinise(*coset_word_acceptor);
    delete coset_word_acceptor;

    State_ID *state = new State_ID[MAX_WORD+1];
    Word_Length depth = 0;
    Node_Count done = 0;
    Node_Count flagged = 0;
    Equation_Word ew(nm);
    int gap = 3;
    ew.reduce();// needed to get a certificate

    state[0] = visible_words->initial_state();
    for (;;)
    {
      Ordinal v = ew.values[depth]+1;
      if (v < pd.coset_symbol)
      {
        ew.set_length(depth);
        ew.append(v);
        State ns = ew.state_word(true);
        State_ID si = state[depth+1] = visible_words->new_state(state[depth],v);
        if (ns && !ns->is_final())
        {
          if (!(char) ++done)
            if (container_status(2,gap,"Looking for coset-inaccessible nodes"
                                 " (" FMT_NC " of " FMT_NC "). "
                                 FMT_NC " found\n",
                                 done,nm.stats.nc[language_L0],flagged))
              gap = 1;
          if (!si)
          {
            ns->node(nm).set_flags(NF_INVISIBLE);
            flagged++;
          }
          depth++;
        }
      }
      else
      {
        if (!depth)
          break;
        depth--;
      }
    }
    delete [] state;
    delete visible_words;
  }
  else if (coset_word_acceptor)
    delete coset_word_acceptor;
}

/**/

bool Rewriter_Machine::schedule_expand(Multi_Node_List & mnl_new,
                                       Multi_Node_List & mnl_old,
                                       Equation_Handle e,bool force)
{
  /* This method decides whether or not to put the specified equation into
     the list of equations that will be expanded on this pass */

  int reason = 0; /* gets set to non-zero value if we decide to expand */
  Word_Length l = e->lhs_length(nm);
  Word_Length r = e->raw_reduced_length(nm);
  Total_Length tl = l+r;

  /* See if one of the criteria for expanding this equation holds.
     The various non-zero values of reason are there just to make
     debugging easier, except that all the values that are != 1
     should stay that way */

  if (tl <= stats.auto_expand_limit)
    reason = 1;
  else if (e->flagged(EQ_AXIOM))
    reason = 2;
  else if (l <= stats.lhs_limit && 2*r <= stats.auto_expand_limit)
    reason = 3;
  else if (tl <= stats.visible_limit && e->flagged(EQ_INTERESTING) && stats.favour_differences)
    reason = 4;
  else if (!pd.inversion_difficult &&
           e->flagged(EQ_HAS_DIFFERENCES+EQ_HAS_PRIMARY_DIFFERENCES) &&
          stats.vital_attained && stats.favour_differences)
    reason = 5;
  else if (maf.options.detect_finite_index==2 && pd.is_coset_system &&
           e->first_letter(nm) >= pd.coset_symbol &&
           e->last_letter() < pd.coset_symbol)
    reason = 6;
  else if (maf.options.strategy & MSF_LEFT_EXPAND &&
          l + 2 <= stats.auto_expand_limit &&
           l >= r && tl <= stats.max_expanded)
    reason = 7;
  else if (force)
    reason = 8;

  if (stats.no_coset_pool &&
      e->raw_reduced_node(nm)->flagged(NF_INVISIBLE) &&
      e->parent(nm)->flagged(NF_INVISIBLE) &&
      e->overlap_suffix(nm)->flagged(NF_INVISIBLE))
    reason = 0;


  if (reason != 0)
  {
    if (reason == 1)
    {
      if (l > stats.lhs_limit)
        stats.lhs_limit = l;
    }
    if (l > stats.max_expanded_lhs)
      stats.max_expanded_lhs = l;
    if (tl > stats.max_expanded)
      stats.max_expanded = tl;
    if (e->status() != NS_Expanded_Equation)
      e->node(nm).set_status(nm,NS_Expanded_Equation);
    e->attach(nm)->set_flags(EQ_EXPANDED);
    if (expand_ont)
      expand_ont->add_state(e);
    else if (e->expanded_timestamp())
      mnl_old.add_state(e);
    else
      mnl_new.add_state(e);
    return true;
  }
  return false;
}

/**/

void Rewriter_Machine::add_to_schedule(Equation_Handle e)
{
  /* if allowable add a newly discovered equation to the schedule instead of
     waiting for the next pass */

  if (maf.options.no_h && e->last_letter() > pd.coset_symbol)
    return;

  if (stats.no_coset_pool &&
      e->raw_reduced_node(nm)->flagged(NF_INVISIBLE) &&
      e->parent(nm)->flagged(NF_INVISIBLE) &&
      e->overlap_suffix(nm)->flagged(NF_INVISIBLE))
    return;

  if (!(maf.options.strategy & MSF_CLOSED) && (expand_ont ? expand_ont->length() : expand_list.length()))
  {
    if (expand_ont)
    {
      e->attach(nm);
      expand_ont->add_state(e);
    }
    else if (!e->flagged(EQ_EXPANDED))
    {
      e->attach(nm);
      expand_list.add_state(e);
    }
    e->node(nm).set_flags(EQ_EXPANDED);
  }
}

/**/

void Rewriter_Machine::coset_extras(Equation_Handle e)
{
  if (pd.is_normal_coset_system && e->is_valid() &&
      e->first_letter(nm) == pd.coset_symbol && !stats.no_coset_pool)
  {
    /* we have an equation that looks like _N*u = _N*v.
       which means that u*V is a member of the subgroup
       Therefore,since _N is normal, we should have _N*g*u*V*G=_N for all g.
       and hence _N*g*u= _N*g*v */
    Word_Length l = e->length(nm);
    Word_Length r = e->reduced_length(nm);
    Ordinal_Word lhs_word(alphabet(),l+1);
    Ordinal_Word rhs_word(alphabet(),r+1);
    lhs_word.set_code(0,pd.coset_symbol);
    rhs_word.set_code(0,pd.coset_symbol);
    e->read_values(lhs_word.data(1),nm,l);
    e->raw_reduced_node(nm)->read_values(rhs_word.data(1),nm);

    for (Ordinal g = 0; g < pd.coset_symbol;g++)
    {
      lhs_word.set_code(1,g);
      rhs_word.set_code(1,g);
      Working_Equation we(nm,lhs_word,rhs_word,Derivation(BDT_Normal_Closure,e,g));
      add_equation(&we,AE_NO_REPEAT);
    }

    update_machine();
  }
}

/**/

void Rewriter_Machine::expand(Node_Handle e1,Node_Count done,unsigned ost)
{
  /* Look for new equations that can be formed by finding
     equations e2 which begin with a suffix of e1.*/
  Node_List candidates_list;
  Simple_Equation se1(nm,e1);
  const Equation_Word lhs_word(nm,e1);
  Ordinal first_letter = se1.first_letter();

  if (maf.options.no_h && e1->last_letter() > pd.coset_symbol)
    return;

  bool first_time = e1->expanded_timestamp()==0;
  Node_Count to_expand = expand_ont ? expand_ont->node_count() :
                         expand_list.length()+re_expand_list.length();
  status(Priority_General,2,1,"Expanding Nodes at LOD %d (" FMT_NC " of "
                              FMT_NC "). Size L/R %d/%d\n",
         ost,done,done+to_expand,se1.lhs_length,se1.rhs_length);
  while (e1->fast_is_primary())
  {
    ID timestamp = e1->expanded_timestamp()+1;

    /* On entry ost will usually be OST_Full, but might be one of the other
       values. This used to be only done when the shallow probe flag was not
       set, but this does not seem sensible */
    if (!maf.options.fast && ost != OST_Aborted_Pass)
    {
      if (e1->flagged(EQ_HAS_DIFFERENCES))
        ost = OST_Full;
      else if (e1->flagged(EQ_INTERESTING) && ost < OST_Moderate)
        ost = OST_Moderate;
      else if (ost == OST_Minimal)
        ost = OST_Quick;
    }
    ID max_id = ost > OST_Minimal ? nm.last_id : (timestamp != 1 ? stats.start_id : 0);
    if (!max_id)
    {
      if (se1.total_length() <= nm.overlap_filter.vital_limit ||
          pd.is_normal_coset_system && first_letter == pd.coset_symbol)
        max_id = nm.last_id;
    }
    if (max_id && maf.options.strategy & MSF_USE_ERAS)
      max_id = stats.start_id;

    if (max_id)
    {
      {
        Find_Candidate find(nm,&candidates_list,timestamp,max_id,
                            e1,maf.options.probe_style & 1 ? stats.max_primary_lhs : stats.max_expanded_lhs,
                            lhs_word,(Overlap_Search_Type) ost,
                            stats.favour_differences);
        Node_Reference node = e1->overlap_suffix(nm);
        node->find_candidates(find,node,se1.lhs_length);
      }
      if (e1->is_adopted())
        e1->node(nm).reduction.expand_timestamp = max_id; /* we set this now in case rhs changes */

      Node_Reference e2;
      int offset;
      Node_Count to_do = candidates_list.length();
      done = 0;
      int gap = 2;

      while (candidates_list.use(&e2,nm,&offset))
      {
        done++;
        if (e1->is_valid())
        {
          Simple_Equation se2(nm,e2);
          consider(e1,e2,Word_Length(offset),se1,se2);
          if (status(Priority_Consider,2,gap,
                     "Considering overlaps (" FMT_NC " of " FMT_NC ") (" FMT_NC " nodes waiting)\n",
                      done,to_do,to_expand))
            gap = 1;
        }
        e2->detach(nm,e2);
      }

      if (e1->is_adopted() && stats.some_ignored)
      {
        e1->node(nm).reduction.expand_timestamp = 0;
        stats.some_ignored = false;
        break;
      }
      if (maf.options.repeat <= 1 || done==0 ||
          e1->expanded_timestamp()==nm.last_id ||
          e1->status() != NS_Expanded_Equation ||
          ost != OST_Full)
        break;
    }
    else
      break;
  }

  if (first_time)
    coset_extras(e1);
}

/**/

int Rewriter_Machine::consider(Equation_Handle e1,Equation_Handle e2,
                               Word_Length offset,
                               const Simple_Equation & eq1,
                               const Simple_Equation & eq2,
                               unsigned flags)
{
  /* Called when expand() has identified that a suffix of e1 and
     a prefix of e2 are the same. This means we can perform two different
     reductions on aXb where aX is the lhs of e1, and Xb is the lhs of e2.
     We form the equation with the "left reduction" cb as the LHS and
     "right reduction" ad as the RHS.
  */
  int retcode = 0;

  if (maf.options.no_h && e1->last_letter() > pd.coset_symbol)
    return 0;

  if (maf.options.max_overlap_length &&
      offset + eq2.lhs_length > maf.options.max_overlap_length)
  {
    stats.some_ignored = true;
    return 0;
  }

  if (e2->is_valid())
  {
    Word_Length common_length = eq1.lhs_length - offset;
    Total_Length l = eq1.rhs_length+eq2.lhs_length-common_length;
    Total_Length r = offset+eq2.rhs_length;

    /* Test equation is not too big to form. We probably must do this now
       even though filtering would probably usually dispose of long equations
       since the RHSes might have changed since the overlap was put on the
       list of overlaps to be considered */
    if (l > MAX_WORD || r > MAX_WORD)
    {
      stats.some_discarded = true;
      return -1;
    }

    bool extended_consider = maf.options.extended_consider && flags == 0;
    if (!flags)
    {
      flags = AE_KEEP_FAILED|AE_UPDATE;
      if (!e2->is_adopted() ||
          maf.options.strategy & MSF_SELECTIVE_PROBE &&
          e2->status() != NS_Expanded_Equation ||
          l+r > stats.visible_limit)
      {
        flags |= AE_DISCARDABLE;
        flags &= ~AE_KEEP_FAILED;
      }
    }

    Ordinal_Word lhs_word(alphabet(),Word_Length(l));
    Ordinal_Word rhs_word(alphabet(),Word_Length(r));
    if (!eq1.rhs_length || !eq2.rhs_length ||
        offset+eq2.lhs_length==eq1.lhs_length)
      flags |= AE_PRE_CANCEL;
    word_copy(lhs_word,eq1.rhs_word,eq1.rhs_length);
    word_copy(lhs_word.data(eq1.rhs_length),
              eq2.lhs_word.data(common_length),
              (eq2.lhs_length-common_length));
    word_copy(rhs_word,eq1.lhs_word,offset);
    word_copy(rhs_word.data(offset),eq2.rhs_word,eq2.rhs_length);

    if (extended_consider)
    {
      Word_DB wdb(alphabet(),1024);
      wdb.insert(lhs_word);
      wdb.insert(rhs_word);
      consider(wdb);
    }
    else
    {
      Working_Equation we(nm,lhs_word,rhs_word,
                          Derivation(BDT_Conjunction,e1,e2,offset));
      retcode = add_equation(&we,flags|AE_NO_REPEAT);
    }
  }
  else
    retcode = -1;
  return retcode;
}

/**/

bool Rewriter_Machine::add_to_reduced_words(Sorted_Word_List &reduced_swl,
                                            Word_DB &next_wdb,
                                            Ordinal_Word * rword)
{
  reduced_swl.insert(*rword);
  if (reduced_swl.count() == 2)
  {
    Working_Equation we(nm,*reduced_swl.word(1),*reduced_swl.word(0),Derivation(BDT_Unspecified));
    bool retcode = add_equation(&we,AE_UPDATE)!=0;
    *rword = *reduced_swl.word(0);
    reduced_swl.empty();
    next_wdb.insert(*rword);
    return retcode;
  }
  return false;
}

/**/

void Rewriter_Machine::consider(Word_DB &current_wdb)
{
  /* This is a method which tries to find more than one reduction of
     each side of an equation that came from KB. There is no reason at all
     why such equations should not themselves contain overlaps, or words
     whose first reduction has an overlap. This method is only used when
     the appropriate tactics flag is set.

     On entry current_wdb contains a list of one or more words that are known
     to be equal as elements.
     We then proceed to try to construct all possible reductions of these words.
     Each time we find we have two different reduced words we form a new
     equation.
     In fact this is not a very viable option, especially for recursive
     orderings, because the list may be astronomically big. So we switch
     to left reduction once we have been going for more than 5 seconds or
     have found more than 2048 words.
  */
  Ordinal_Word test_word(alphabet());
  Ordinal_Word reduced_word(alphabet());
  Word_DB next_wdb(alphabet(),1024);
  Word_DB total_wdb(alphabet(),1024);
  Sorted_Word_List reduced_swl(alphabet());
  Element_Count count;
  Progress_Count start = stats.status_count;
  bool extended = true;
  int gap = 2;

  while ((count = current_wdb.count())!=0)
  {
    for (Element_ID word_nr = 0; word_nr < count;word_nr++)
    {
      current_wdb.get(&test_word,word_nr);
      Word_Length l = test_word.length();
      Ordinal * values = test_word.buffer();
      Node_Reference s = nm.start();
      bool found = false;
      /* Find all possible first reductions of the current entry */
      for (Word_Length valid_length = 0; valid_length < l;)
      {
        s = s->transition(nm,values[valid_length++]);
        if (s->is_equation())
        {
          s = s->primary(nm,s);
          const Node * rs = s->fast_reduced_node(nm);
          Word_Length lhs_length = s->lhs_length(nm);
          Word_Length rhs_length = rs->length(nm);
          found = true;
          if (l-lhs_length + rhs_length <= MAX_WORD)
          {
            reduced_word.allocate(l - lhs_length + rhs_length,false);
            Ordinal * new_values = reduced_word.buffer();
            word_copy(new_values,values,valid_length-lhs_length);
            rs->read_values(new_values+valid_length-lhs_length,nm,rhs_length);
            word_copy(new_values+valid_length-lhs_length+rhs_length,
                      values+valid_length,l-valid_length);

            if (extended)
            {
              if (status(Priority_Extended,2,gap,
                         "Extended consideration of overlap\n"))
                gap = 1;
            }
            else if (status(Priority_Extended,2,gap,
                            "Normal consideration of overlap\n"))
              gap = 1;
            if (lhs_length == l && !rs->is_final())
              add_to_reduced_words(reduced_swl,next_wdb,&reduced_word);
            else
            {
              if (!total_wdb.contains(reduced_word))
              {
                next_wdb.insert(reduced_word);
                if (extended)
                  total_wdb.insert(reduced_word);
              }
              valid_length -= lhs_length-1;
              s = nm.start();
              if (count > 2048 || stats.status_count-start > 5)
              {
                extended = false;
                break;
              }
              if (l > stats.visible_limit/2 || !extended)
                break;
            }
          }
          else
          {
            valid_length -= lhs_length-1;
            s = nm.start();
          }
        }
      }
      if (!found)
        add_to_reduced_words(reduced_swl,next_wdb,&test_word);
    }
    current_wdb.take(next_wdb);
  }
}

/**/

bool Rewriter_Machine::add_axiom(const Word & lhs,const Word & rhs,unsigned flags)
{
  /* As well as proper axioms, this method is used to set up equations
     from pkbprog and kbprog files when the -r option is used.
  */
  if (pd.is_coset_system)
    generator_properties[pd.coset_symbol] |= GF_RIGHT_CANCELS;
  Working_Equation we(nm,lhs,rhs,Derivation(flags & AA_RESUME ? BDT_Known : BDT_Axiom));

  /* If none of the axioms have odd length then no equation can have
     odd length, so limits should go up in steps of 2 */
  if (we.total_length() & 1)
    stats.total_increment = 1;

  if (flags & AA_DEDUCE)
  {
    if (!stats.started)
      start();
    if (flags & AA_ELIMINATE)
    {
      Total_Length ael = stats.auto_expand_limit;
      Total_Length vl = stats.visible_limit;
      if (!we.normalise(AE_SAFE_REDUCE|AE_REDUCE))
        return false;
      stats.auto_expand_limit = we.total_length()+4;
      stats.visible_limit = stats.auto_expand_limit+8;
      int i = 0;
      int j = nr_generators;
      while (j > 4 && i < 4)
      {
        stats.auto_expand_limit--;
        stats.visible_limit--;
        j = (j+1)/2;
        i++;
      }
      if (stats.auto_expand_limit < ael)
        stats.auto_expand_limit = ael;
      if (stats.visible_limit < vl)
        stats.visible_limit = vl;
      if (stats.visible_limit > 20)
      {
        stats.visible_limit = stats.auto_expand_limit = we.total_length();
        if (stats.visible_limit < 20)
          stats.visible_limit = 20;
      }
      expand_machine(true);
    }
  }
  unsigned ae_flags = flags & AA_DEDUCE ? AE_UPDATE : 0;
  ae_flags |= AE_NO_REPEAT;
  if (flags & AA_POLISH)
    ae_flags |= AE_SAVE_AXIOM;
  if (flags & AA_RESUME)
  {
    if (!stats.want_differences)
      ae_flags |= AE_INSERT;
    else
      ae_flags |= AE_KEEP_LHS;
    stats.resumed = true;
  }
  bool save_throttle = maf.options.no_throttle;
  if (flags & AA_RESUME)
    maf.options.no_throttle = true;
  if (add_equation(&we,ae_flags)==1)
  {
    maf.options.no_throttle = save_throttle;
    if (flags & AA_DEDUCE)
    {
      set_initial_special_limit();
      update_machine(UM_CHECK_POOL);
    }
    /* I also used to call update_machine() after inserting each equation from
       a resumed run. But this may mean it takes a very long time even to
       get started. So now there is just one massive update_machine()
       when we we start KB */

    return true;
  }
  maf.options.no_throttle = save_throttle;
  if ((flags & AA_WARN_IF_REDUNDANT+AA_RESUME) == AA_WARN_IF_REDUNDANT)
  {
    if (container.progress(1,"Redundant axiom:"))
    {
      we.lhs_word() = lhs;
      we.rhs_word() = rhs;
      we.print(container,container.get_log_stream());
    }
  }
  return false;
}

/**/

int Rewriter_Machine::add_equation(Working_Equation * we,unsigned flags)
{
  /* This method tries to add an equation ax=b to the tree.
     The equation may well instead be pooled, or even discarded altogether
     if Equation_Manager does not like the look of it.
     Even if insertion succeeds, afterwards ax and b may not be regarded as
     equal, because the equation might have been balanced past an equation
     bX=c for which cx=b had not been created (normally only possible with
     recursive orderings).
     So this method optionally keeps inserting the same equation, until
     it actually is recognised that ax=b, or learn() does not return 1.
  */

  int retcode = 0;
  if (flags & AE_NO_REPEAT)
    retcode = em.learn(we,flags,&stats.some_discarded);
  else
  {
    we->failed = true;
    we->changed = false;
    for (;;)
    {
      Working_Equation we2(*we);
      int local_retcode = em.learn(&we2,flags,&stats.some_discarded);
      if (!we2.failed)
        we->failed = false;
      if (retcode==0 || local_retcode==1)
        retcode = local_retcode;
      if (local_retcode != 0)
        if (we2.changed)
          we->changed = true;
      if (local_retcode != 1 || !we2.balanced)
        break;
      flags |= AE_NO_UNBALANCE;
    }
  }

  if (retcode==1 && flags & AE_UPDATE)
  {
    // now optionally find its "obvious" consequences

    /* we first check the pool is not getting too big in comparison with
       the tree. This used to happen a lot, but now that most of the
       equations that would have gone into the pool are thrown away it
       is very unlikely */
    if (pool && pool->count() > 200000 &&
        jm.status().nr_visible < (Node_Count) pool->count()/20 &&
        !stats.pool_full)
    {
      if (nr_good_generators > 1)
      {
        Element_Count count = pool->count();

        while (count/nr_good_generators > 100)
        {
          count /= nr_good_generators;
          stats.visible_limit += stats.total_increment;
        }
      }
      else
        stats.visible_limit += 2;
      if (stats.visible_limit > stats.auto_expand_limit*2-2)
      {
        stats.visible_limit = stats.auto_expand_limit*2-2;
        stats.pool_full = true;
      }
      flags |= UM_CHECK_POOL;
    }
    update_machine(flags);
  }
  return retcode;
}

/**/

void Rewriter_Machine::add_to_tree(Working_Equation * we,unsigned flags,
                                   Total_Length was_total)
{
  /* On entry we is assumed to be a properly insertable equation ug=v:
     i.e. the lhs has an irreducible prefix, the rhs is irreducible and
     lhs > rhs.

     The equation might not be balanced, i.e. it might be the case that
     u > vG.

     It is absolutely vital that we call Node_Manager::check_transitions()
     if we create a new LHS or make a previously irreducible word into a
     prefix, and that we don't try to use the tree to perform reductions
     in between doing the former and the latter.

     When this method returns it is guaranteed that there are no links
     into removed nodes, except possibly for the RHS of some equations.
  */

  /* The trace has to have the equation before the tree changes, since
     what it cites is resolved against the tree */
  int pid = proof_trace ? proof_trace->emit_equation(*we) : 0;
  Equation_Word & lhs = we->lhs_word();
  Equation_Word & rhs = we->rhs_word();
  Word_Length l = lhs.word_length;
  Word_Length r = rhs.word_length;
  Total_Length tl = l+r;
  Node_Reference old_e = lhs.state_lhs(false);
  Node_Reference e;
  Node_Reference prefix = lhs.state_prefix();
  Node_Reference rhs_node = rhs.state_word(true);
  Node_Reference trailing_subword;
  Word_Length new_prefix_length = l;
  bool must_simplify = (flags & AE_MUST_SIMPLIFY)!=0;
  Ordinal lvalue = lhs.first_letter();
  Ordinal rvalue = lhs.last_letter();
  bool new_reduction = true;

  if (old_e->length(nm) == lhs.word_length)
    e = old_e;

  trailing_subword = old_e->primary(nm,old_e);
  if (trailing_subword == e)
    if (e->is_final())
      trailing_subword = e->overlap_suffix(nm);
    else
      trailing_subword = Node_Reference(0,0);

  if (rhs_node.is_null())
    rhs_node = rhs.get_node(r,false,0);
  if (old_e->is_final())
    new_reduction = false;
  if (new_reduction || maf.options.secondary>=3)
    must_simplify = true;

  if (prefix.is_null() || must_simplify && !prefix->is_prefix())
    prefix = lhs.get_node(l-1,must_simplify,
                          must_simplify ? &new_prefix_length : 0);

  if (pd.is_coset_system && lvalue <= pd.coset_symbol && new_reduction &&
      !prefix->flagged(NF_INVISIBLE))
    stats.visibility_correct = false;

  if (trailing_subword.is_null())
    trailing_subword = nm.get_state(lhs,l,1);
  if (!e.is_null() && e->is_final())
    jm.cancel_jobs(e);

  if (r==0)
  {
    generator_properties[lvalue] |= GF_RIGHT_CANCELS|GF_RIGHT_INVERSE;
    generator_properties[rvalue] |= GF_LEFT_CANCELS|GF_LEFT_INVERSE;
    if (l == 1)
    {
      /* the if below should be unnecessary, but we want to make absolutely
         sure we don't miscount the number of trivial generators, as
         reduce() will use this to optimise reduction in the case all
         generators become trivial */
      if (!(generator_properties[lvalue] & GF_TRIVIAL))
      {
        nr_trivial_generators++;
        generator_properties[lvalue] |= GF_TRIVIAL;
      }
    }
  }

  if (l == 1)
  {
    if (!redundant(rvalue))
    {
      nr_good_generators--;
      generator_properties[rvalue] |= GF_REDUNDANT;
      container.progress(1,"Generator %s is being eliminated!\n",
                         alphabet().glyph(rvalue).string());
    }

  }
  if (pd.is_coset_system)
  {
    if (l == 2 && lvalue == pd.coset_symbol)
      if (!(generator_properties[rvalue] & GF_COSET_REDUCIBLE))
      {
        nr_coset_reducible_generators++;
        generator_properties[rvalue] |= GF_COSET_REDUCIBLE;
      }
  }

  if (we->eq_flags & EQ_AXIOM)
  {
    if (!stats.started || flags & AE_SAVE_AXIOM)
    {
      Total_Length safe = min(tl,was_total);
      if (tl > stats.max_expanded)
      {
        stats.max_expanded = tl;
        if (safe > stats.auto_expand_limit)
          if (!(pd.is_coset_system && lvalue >= pd.coset_symbol))
            stats.auto_expand_limit = safe;
          else
            safe = 0;
      }
      if (safe+4 > stats.visible_limit)
        stats.visible_limit = safe+4;
      if (l > stats.lhs_limit)
        stats.lhs_limit = l;
    }
  }
  else if (maf.options.log_level >=2)
  {
    bool pure = true;
    for (Word_Length i = 1;i < l;i++)
      if (we->lhs->values[i] != we->lhs->values[0])
      {
        pure = false;
        break;
      }

    if (pure || r == 0 || pd.is_coset_system && l==2 && lhs.value(0)==pd.coset_symbol)
      we->print(container,container.get_log_stream());
  }

  if (!e.is_null())
  {
    if (!e->is_final())
    {
      if (e->flagged(NF_IS_DIFFERENCE+NF_IS_HALF_DIFFERENCE) &&
          !rhs_node->reduced.inverse && !maf.options.no_half_differences)
        rhs_node->node(nm).set_inverse(nm,rhs_node);
      Ordinal rhs_last;
      if (e->lhs_node && right_cancels(rhs_last = rhs_node->last_letter()))
      {
        Node_Reference other;
        for (other = e->first_lhs(nm);other;
             other = other->next_lhs(nm))
          if (rhs_last == other->last_letter())
          {
            /* from ax=by and cy=ax deduce c=b */
            Derivation d(BDT_Equal_RHS,other,e);
            if (proof_trace)
              d.set_trace(proof_trace->node_equation(nm,other),pid);
            jm.equate_nodes(other->parent(nm),rhs_node->parent(nm),d);
            break;
          }
      }
    }
    else
    {
      if (e->flagged(EQ_INTERESTING+EQ_HAS_DIFFERENCES+EQ_HAS_PRIMARY_DIFFERENCES))
        remove_differences(e);
      if (flags & AE_STRONG_INSERT)
      {
        Derivation d(BDT_Equal_RHS,e,e);
        if (proof_trace)
          d.set_trace(proof_trace->node_equation(nm,e),pid);
        jm.equate_nodes(e->fast_reduced_node(nm),rhs_node,d);
        flags &= ~AE_STRONG_INSERT;
      }
    }
  }
  e = nm.construct_equation(prefix,rvalue,we->eq_flags,rhs_node,trailing_subword);
  lhs.state[l] = e;
  if (proof_trace)
    proof_trace->set_node(e,pid);

  /* Having constructed ax=b we now think about bX=a, which should
     be an equation as well.*/
  bool urgent = true;
  bool balanced = true;
  Ordinal ivalue = inverse(rvalue);
  if (ivalue != INVALID_SYMBOL)
  {
    Node_Reference conjugate = rhs_node->transition(nm,ivalue);
    if (conjugate->reduced_node(nm,conjugate) != prefix && r < MAX_WORD)
    {
      if (!pd.g_level)
      {
        rhs.append(ivalue);
        urgent = balanced = rhs.compare(Subword(lhs,0,l-1)) > 0;
        rhs.set_length(r);
      }
      Derivation d(BDT_Right_Conjugate,e);
      if (proof_trace)
        d.set_trace(pid,0);
      jm.equate_nodes(rhs_node,prefix,d,
                      ivalue,urgent && !(flags & AE_WANT_WEAK_INSERT));
    }
  }

  const Language language = e->language(nm);
  if (language==language_L1)
  {
    if (tl < stats.shortest_recent && tl > 4)
      stats.shortest_recent = tl;
    stats.equation_delta++;
    stats.want_overlaps = true;
    stats.complete = stats.primary_complete = false;
    must_simplify = true;
    if (pd.is_coset_system)
    {
      if (lvalue < pd.coset_symbol)
        stats.g_complete = false;
      else if (lvalue > pd.coset_symbol && rvalue > pd.coset_symbol)
        stats.h_complete = false;
      else if (new_reduction)
      {
        stats.coset_complete = false;
        stats.last_coset_equation = stats.status_count;
      }
    }
  }

  if (maf.options.log_flags & LOG_EQUATIONS)
  {
    /* We have to do this before the equations involved in the
       deduction have been changed */
    we->print_derivation(container.get_log_stream());
  }
  if (maf.options.special_overlaps == 5)
  {
    new_prefix_length = 0;
    must_simplify = true;
  }
  if (must_simplify) /* Ensure no links to removed nodes survive */
    nm.check_transitions(lhs,new_prefix_length,true,stats.want_overlaps,new_reduction);

  /* Look for deductions */
  if (e->lhs_node)
  {
    if (right_cancels(rvalue))
    {
      for (Node_Reference other = e->next_lhs(nm);other;other = other->next_lhs(nm))
        if (other->last_letter() == rvalue)
        {
          /* from ax=b and cx=b deduce a=c */
          Derivation d(BDT_Equal_RHS,other,e);
          if (proof_trace)
            d.set_trace(proof_trace->node_equation(nm,other),pid);
          jm.equate_nodes(other->parent(nm),prefix,d);
          break;
        }
    }
    if (left_cancels(lvalue) && language==language_L1)
    {
      for (Node_Reference other = e->next_lhs(nm);other;other = other->next_lhs(nm))
        if (lvalue == other->first_letter(nm) && other->fast_is_primary())
        {
          /* from ax=b a==gu and gv=b deduce ux=v */
          Derivation d(BDT_Equal_RHS,other,e);
          if (proof_trace)
            d.set_trace(proof_trace->node_equation(nm,other),pid);
          jm.equate_nodes(other,e,d);
          break;
        }
    }
  }

  jm.schedule_jobs(e,choose_status(e,!balanced ||
                                   stats.started && maf.options.strategy & MSF_CONJUGATE_FIRST ?
                                   NS_Unconjugated_Equation :
                                   NS_Undifferenced_Equation));
  if (e != old_e && !e->fast_is_primary() && !old_e->fast_is_primary())
    jm.schedule_overlap(e,old_e,l - old_e->length(nm),true);
#ifdef DEBUG
  const Equation_Stats & estats = jm.status();
  if (estats.nr_visible != nm.stats.ntc[1] ||
      nm.stats.nc[0] != nm.stats.ntc[0]+1 ||
      estats.nr_visible != nm.stats.nc[1]+nm.stats.nc[2]+nm.stats.nc[3])
    MAF_INTERNAL_ERROR(container,("Node_Manager and Job_Manager inconsistent!\n"));
#endif
  if ((!stats.want_differences || maf.options.is_kbprog) && stats.started &&
      (maf.options.max_equations && nm.stats.nc[1] >= maf.options.max_equations ||
       maf.options.max_time && stats.status_count > maf.options.max_time))
    maf.aborting = true;
}

/**/

bool Rewriter_Machine::add_to_pool(Element_ID * id,Working_Equation *we,unsigned flags)
{
  bool good = false;
  Equation_Word & lhs = we->lhs_word();
  Equation_Word & rhs = we->rhs_word();
  Word_Length l = lhs.word_length;
  Word_Length r = rhs.word_length;
  Total_Length tl = l+r;

  if (we->failed &&
      tl > stats.visible_limit && !(flags & AE_KEEP_FAILED))
    return false;

  if (!pool)
  {
    unsigned suggested_size = jm.status().nr_equations - jm.status().nr_visible;
    if (suggested_size < 1024)
      suggested_size = 1024;
    pool = new Equation_DB(alphabet(),suggested_size);
  }

  if (lhs.language == language_L0)
  {
    stats.complete = false;
    good = !we->failed;
    if (!we->failed)
    {
      if (tl < stats.pool_limit)
      {
        stats.pool_limit = tl;
        stats.discard_limit = stats.pool_limit + POOL_DEPTH;
      }
      if (l <= stats.best_pool_lhs)
      {
        if (l < stats.best_pool_lhs || tl < stats.best_pool_lhs_total)
          stats.best_pool_lhs_total = tl;
        stats.best_pool_lhs = l;
      }
      if (r <= stats.best_pool_rhs)
      {
        if (r < stats.best_pool_rhs || tl < stats.best_pool_rhs_total)
          stats.best_pool_rhs_total = tl;
        stats.best_pool_rhs = r;
      }
    }

    if (pd.is_coset_system)
    {
      Ordinal lvalue = lhs.first_letter();
      if (lvalue < pd.coset_symbol)
        stats.g_complete = false;
      else if (lvalue > pd.coset_symbol)
        stats.h_complete = false;
      else
        stats.coset_complete = false;
    }
  }

  if (dt && !stats.favour_differences && !stats.hidden_differences)
    stats.hidden_differences = dt->is_interesting_equation(*we,AE_IGNORE_PRIMARY);

  /* Pool equations are found again by their words, see maf_trace.cpp */
  int pid = proof_trace ? proof_trace->emit_equation(*we) : 0;
  Element_ID key;
  if (!pool->insert(lhs,&key))
  {
    *id = key;
    Ordinal_Word old_rhs_word(alphabet());
    pool->get_rhs(&old_rhs_word,key);
    if (words_differ(rhs,old_rhs_word,r+1))
    {
      Derivation d(key,BDT_Equal_LHS);
      if (proof_trace)
        d.set_trace(proof_trace->equation(lhs,old_rhs_word),pid);
      em.queue(old_rhs_word,rhs,d,AE_PRE_CANCEL);
      if (old_rhs_word.compare(rhs) > 0)
      {
        if (good && !stats.want_differences)
          stats.equation_delta++;
        if (maf.options.log_flags & LOG_EQUATIONS)
          we->print_derivation(container.get_log_stream(),key);
        pool->update_rhs(key,rhs);
        return true;
      }
    }
    return false; /* The same or a better equation is already in the pool */
  }
  else
  {
    *id = key;
    jm.change_pool_count(1);
    if (good && !stats.want_differences)
      stats.equation_delta++;
    if (maf.options.log_flags & LOG_EQUATIONS)
      we->print_derivation(container.get_log_stream(),key);
    pool->update_rhs(key,rhs);
  }
  if (tl < stats.shortest_recent)
    stats.shortest_recent = tl;
  return true;
}

/**/

void Rewriter_Machine::schedule_optimise(const Word &word)
{
  jm.schedule_optimise(word);
}

/**/

void Rewriter_Machine::optimise(const Word & word)
{
  Word_Length l = word.length();
  Equation_Word lhs_word(nm,word);
  lhs_word.set_length(0);
  const Ordinal * values = word.buffer();
  int gap = 2;
  for (Word_Length i = 0; i < l;i++)
  {
    if (status(-1,2,gap,"Optimising slow reduction (%d of %d steps)\n",i,l))
      gap = 1;
    lhs_word.append(values[i]);
    if (lhs_word.language != language_L0)
    {
      Node_Reference e = lhs_word.state_lhs(true);
      if (e.is_null())
      {
        bool failed = false;
        lhs_word.reduce(AE_KEEP_LHS,&failed);
        Working_Equation we(nm,lhs_word);
        if (add_equation(&we,AE_KEEP_LHS|AE_COPY_LHS|AE_NO_UNBALANCE|AE_NO_REPEAT|AE_NO_BALANCE|AE_WANT_WEAK_INSERT) == 1)
        {
          if (lhs_word.compare(we.lhs_word()) == 0)
            lhs_word = we.rhs_word();
        }
        else if (we.lhs_word().fast_length() &&
                 we.lhs_word().state_lhs(true).is_null())
          break;
        if (we.failed)
          break;
      }
      else
        lhs_word.read_word(e->fast_reduced_node(nm),true);
    }
  }
}

/**/

void Rewriter_Machine::set_initial_special_limit()
{
  stats.initial_special_limit = stats.visible_limit;
  if (stats.initial_special_limit > 2*stats.auto_expand_limit)
    stats.initial_special_limit = 2*stats.auto_expand_limit;
  /* We try to keep special overlaps genuinely special, so that they don't
     dominate processing. For wreath product word orderings we allow more of
     them because this seems beneficial, but we restrict anyway if
     MSF_NO_FAVOUR_SHORT is set in any case, because then visible_limit is very
     high */

  if ((!pd.inversion_difficult || maf.options.strategy & MSF_NO_FAVOUR_SHORT) &&
      stats.auto_expand_limit < stats.initial_special_limit)
  {
    if (maf.options.special_overlaps == 1)
    {
      Element_Count i = nr_good_generators;
      while (i > 3 && stats.initial_special_limit >= stats.auto_expand_limit)
      {
        i /= 2;
        stats.initial_special_limit = (stats.auto_expand_limit + 3*stats.initial_special_limit)/4;
      }
    }
    if (stats.total_increment == 2 && stats.initial_special_limit & 1)
      stats.initial_special_limit++;
  }
  if (stats.initial_special_limit > 100)
    stats.initial_special_limit = max(stats.auto_expand_limit,Total_Length(100));
  if (stats.initial_special_limit > height()*3)
    stats.initial_special_limit = height()*3; // in case no_pool is in effect
}

/**/

void Rewriter_Machine::update_machine(unsigned flags)
{
  stats.special_limit = stats.initial_special_limit;
  stats.update_time = stats.status_count;
  jm.update_machine(flags);
}

/**/

bool Rewriter_Machine::build_allowed(bool want_build)
{
  bool retcode = true;
  bool log = want_build;

  if (stats.complete || maf.options.assume_confluent ||
      stats.primary_complete && !pool)
    retcode = false;

  if (retcode && dt && !maf.options.fast)
  {
    bool again = want_build;
    for (;;)
    {
      const Difference_Tracker::Status & dstats = dt->status();
      Total_Length vital = dt->vital_eqn_limit();
      /* This test is also heuristic, but works fairly well. We want to
         be sure that we can at least correctly reduce all words that might
         reduce to word-differences. And we want the size of the equations
         we might generate to have attained the size of the longest difference
         seen so far. However we relax this, and only insist that
         one or other condition probably holds. */
      if ((stats.max_expanded <= dstats.max_difference_eqn &&
           stats.max_expanded <= dstats.max_primary_difference_eqn ||
           nm.overlap_filter.forbidden_limit*2-4 <=
           dstats.max_difference_eqn) && stats.visible_limit < vital)
      {
        if (again)
        {
          explore_dm(false);
          check_differences(CD_CHECK_INTEREST);
          vital = dt->vital_eqn_limit();
          again = false;
          continue;
        }
        if (log)
          container_status(1,1,
                           "Further differences expected based on equation sizes.\n");
        retcode = false;
      }
      break;
    }
  }
  return retcode;
}

/**/

bool Rewriter_Machine::examine(unsigned rm_state)
{
  /* After a pass of Knuth Bendix decide whether to carry on or to
     try to build the automata.
     The return code is true if it might be a good idea to build the automata
     now, and false if KB expansion should continue */

  is_dm_changing(true,false,rm_state==0);
  if (maf.options.fast && rm_state==0 &&
      stats.status_count - stats.last_difference >= stats.normal_timeout*2)
    stats.difference_delta = 0;
  bool retcode = !maf.options.assume_confluent && stats.want_differences &&
                 (rm_state==0 && !stats.difference_delta && pool && stats.equation_delta > 10 ||
                  rm_state & RMS_Hurrying);
  if (rm_state != RMS_Normal)
    stats.complete = stats.approved_complete = false;
  if (rm_state & RMS_Aborting)
    return true;
  if (stats.want_differences)
  {
    if (maf.options.min_time && stats.status_count < maf.options.min_time)
      return false;
    if (maf.options.max_equations && nm.stats.nc[1] > maf.options.max_equations)
      return true;
    if (maf.options.max_time && stats.status_count > maf.options.max_time)
      return true;
  }

  if (rm_state & RMS_Too_Tight+RMS_Aborted_Pass)
    return false;

  if (rm_state == 0 && stats.want_differences)
  {
    if (!stats.difference_delta && !stats.complete && stats.equation_delta)
    {
      if (stats.equation_delta && !maf.options.assume_confluent && pool &&
          stats.upper_consider_limit > dt->status().max_difference_eqn)
        retcode = true;
      if (stats.visible_limit > dt->status().max_interesting_eqn)
        retcode = true;
    }
    if (!stats.favour_differences && stats.hidden_differences)
      retcode = false;
  }

  if (retcode)
    retcode = build_allowed(true);
  return retcode;
}

/**/

void Rewriter_Machine::purge()
{
  /* Get rid of any spare nodes and node list entries. There may be a lot
     of these and this may free up a substantial amount of memory */
  Node_List::purge();
}

/**/

void Rewriter_Machine::kill_node(Node * node)
{
  // This is the last and final call for this node, which is absolutely
  // dead and gone and likely to be re-used.
  if (derivation_db && node->flagged(NF_REDUCIBLE))
  {
    // We don't attach to equations when we stick them in this database,
    // and there is no obvious time when we can delete them, since a
    // deduction might not be processed until an equation that was used
    // in its derivation had been eliminated. So we are going to live with
    // some equations no longer being available
    Element_ID id = derivation_db->find_entry(&node,sizeof(State),false);
    derivation_db->remove_entry(id);
  }
}

/**/

bool Rewriter_Machine::add_correction(const Ordinal_Word & lhs,
                                       const Ordinal_Word & rhs,bool is_primary)
{
  Working_Equation we(nm,lhs,rhs,Derivation(BDT_Diff_Reduction));

  unsigned ae_flags = AE_KEEP_LHS|AE_DISCARDABLE|AE_CORRECTION|AE_NO_REPEAT;
  if (!is_primary)
    ae_flags |= AE_IGNORE_PRIMARY;
  if (pd.is_coset_system ? pd.is_coset_finite : pd.is_g_finite)
    ae_flags |= AE_INSERT;

  if (add_equation(&we,ae_flags)==1)
  {
    if (container.progress(3,"Correcting equation for automata\n"))
    {
      we.print(container,container.get_log_stream());
      container.progress(3,"\n");
    }
    update_machine(0);
    return true;
  }
  return false;
}

/**/

void Rewriter_Machine::clear_differences(Node_List * equation_list)
{
  /* Delete the old difference tracker and prepare to rebuild the new one
     This code used to clear the difference flags as it scanned the nodes
     to build up the equation list. But because the code now ignores
     NF_INVISIBLE nodes this is not correct since these might have had
     differences. The Difference_Tracker will now clear the flags that
     are set already. This will cause flags to be cleared multiple times,
     but since this is very fast it should not matter.
  */
  update_machine(0); // we must make sure there is no outstanding work.
                     // otherwise word-differences that are scheduled for
                     // correction will cause the new tracker to detach
                     // from things it is not attached to.
  hide_equations(0);
  State_Count count = dt->dt_stats.nr_differences;
  delete dt;
  dt = Difference_Tracker::create(nm,count);

  /* initialise the list of equations to process */
  const Equation_Stats &estats = jm.status();
  Multi_Node_List mnl(nm,height());
  {
    Node_Reference nr;
    Node_Count done = 0;
    bool inside = true;
    Node_Iterator ni(nm);
    while (ni.scan(&nr,inside))
    {
      Node & node = nr->node(nm);
      if (node.is_final())
      {
        /* Formerly this code did not attach the equations it put on the
           list it builds on the grounds that the equations could not
           disappear during the lifetime of the list. But this is no longer
           always the case so we now do */
        if (equation_list &&
            (!pd.is_coset_system || node.last_letter() < pd.coset_symbol))
        {
          node.attach(nm);
          mnl.add_state(nr);
        }
        status(Priority_General,2,1,"Cleaning word-differences. (" FMT_NC " of " FMT_NC ")\n",
               ++done,estats.nr_visible);
      }
      else
      {
        if (pd.is_coset_system)
          inside = !node.flagged(NF_INVISIBLE) && ni.word()[0] <= pd.coset_symbol;
      }
    }
  }
  if (equation_list)
    mnl.collapse(equation_list);
}

bool Rewriter_Machine::check_differences(unsigned flags)
{
  if (!stats.want_differences)
    return false;

  if (flags & CD_STABILISE &&
      !(stats.complete || pd.is_coset_system && stats.coset_complete))
  {
    for (;;)
    {
      if (!dt->is_valid(false) || (flags & CD_EXPLORE))
        explore_dm(true);
      dt->scan_all_differences();
      if (!jm.work_pending())
        break;
      update_machine(0);
      if (maf.aborting)
        break;
    }
  }

  update_machine(UM_CHECK_OVERSIZED);

  if (flags & CD_FORCE_REBUILD ||
      !dt->is_valid((flags & CD_CHECK_INTEREST)!=0) ||
      !dt->is_inverse_complete())
  {
    for (;;)
    {
      Node_List equation_list;
      clear_differences(&equation_list);
      /* Insert a state for every equation's word-differences */
      {
        Node_Reference e;
        Node_Count to_do = equation_list.length();
        Node_Count done = 0;
        while (equation_list.use(&e,nm))
        {
          done++;
          if (e->is_valid())
          {
            Working_Equation we(nm,e,e->fast_reduced_node(nm),Derivation(BDT_Known,e));
            dt->learn_equation(we,e,0);
            if (flags & CD_UPDATE)
              update_machine(0);
            container_status(2,1,"Examining equations for word-differences. ("
                                 FMT_NC " of " FMT_NC ")\n"
                                 FMT_ID "/" FMT_ID
                                 " primary/secondary differences so far.\n",
                             done,to_do,
                             dt->dt_stats.nr_primary_differences,
                             dt->dt_stats.nr_differences);
          }
          e->detach(nm,e);
        }
      }

      dt->close();
      dt->take_changes();

      if (jm.work_pending())
      {
        /* We found more equations while we rebuilt the tracker. We will live
           with this unless we want to build the automata now */
        Certificate c = nm.current_certificate();
        stats.difference_delta++;
        update_machine(0);
        if (nm.check_certificate(c) != CV_VALID &&
            flags & CD_STABILISE && !maf.aborting &&
            !dt->is_valid((flags & CD_CHECK_INTEREST)!=0))
          continue;
      }
      break;
    }

    Difference_Tracker::Status dstats = dt->status();
    container.status(1,1,FMT_ID " differences, (" FMT_ID " primary), "
                     "max eqn size=%d,max length=%d.\n",
                      dstats.nr_differences,dstats.nr_primary_differences,
                      dstats.max_primary_difference_eqn,
                      dstats.max_difference_length);
  }
  else
    dt->clear_changes(false);
  return true;
}

/**/

FSA_Simple * Rewriter_Machine::grow_wd(unsigned cd_flags,unsigned gwd_flags)
{
  if (!check_differences(cd_flags))
    return 0;
  return dt->grow_wd(gwd_flags);
}

/**/

bool Rewriter_Machine::explore_acceptor(const FSA *wa,
                                         const FSA * dm2,
                                         bool finite,
                                         Progress_Count allow)
{
  /* This method uses an equation_word() to explore all the nodes in the
     Rewriter_Machine, looking for nodes that correspond to words that
     are rejected by the word-acceptor. If such a node is found the
     missing equation is created to remove the word.
     This method is more or less superseded by explore_dm(), but could
     be useful for largish finite groups with a lot of secondary differences
  */
  State_ID *state = new State_ID[MAX_WORD+1];
  Word_Length depth = 0;
  Word_Length limit = 0;
  bool found = false;
  bool finished = false;
  bool retcode = false;
  bool done = false;
  Progress_Count started = stats.status_count;
  Equation_Word ew(nm);
  Diff_Reduce dr(dm2);
  Working_Equation we(nm,Derivation(BDT_Diff_Reduction));

  state[0] = wa->initial_state();
  if (!allow)
    allow = (finite ? 40 : 20);
  ew.set_length(0);
  ew.reduce();
  for (;;)
  {
    Equation_Word & lhs = we.lhs_word();
    Equation_Word & rhs = we.rhs_word();
    Ordinal v = ew.values[depth]+1;
    if (v < nr_generators)
    {
      ew.set_length(depth);
      ew.append(v);
      State ns = ew.state_lhs(false);
      if (!ns)
      {
        /* As we navigate through the Rewriter_Machine, our equation_word
           will normally have a state. However, if we have inserted an
           equation, this may have become invalid, so we need to call
           reduce() to get the state back. The word should never actually be
           reducible, if it is, then we have discovered so many new equations
           that our acceptor is now essentially worthless as the Rewriter_Machine
           is rejecting words that it accepts */
        if (ew.reduce(AE_KEEP_LHS))
        {
          /* give up - our acceptor has accepted a now reducible word */
          delete [] state;
          return true;
        }
        ns = ew.state_lhs(false);
      }
      if (depth == 0 || v != inverse(ew.values[depth-1]))
      {
        State_ID wa_si = state[depth+1] = wa->new_state(state[depth],v);
        if (!wa_si)
        {
          if (ns &&
              Node::fast_find_node(nm,ew.state[depth+1-limit])->length(nm)+limit == depth+1)
          {
            rhs = lhs = ew;
            bool wrong = !rhs.reduce();
            if (!stats.complete)
            {
              Ordinal_Word rhs_word(rhs);
              wrong = dr.reduce(&rhs_word,rhs_word,0,wa) != 0;
              if (wrong)
              {
                rhs = rhs_word;
                rhs.invalidate();
              }
            }
            if (add_equation(&we,AE_KEEP_LHS|AE_DISCARDABLE|AE_IGNORE_PRIMARY|AE_NO_REPEAT)==1)
              done = true;
          }
          if (status(Priority_General,2,1,"Examining word-acceptor at depth %d pass %d\n",depth,limit))
          {
            if (done)
            {
              update_machine();
              done = false;
            }
            Progress_Count used = stats.status_count - started;
            if (!found && used > allow || used/2 > allow)
              break;
          }
        }
        else if (ns && ns->length(nm) == depth+1 && !ns->is_final())
          depth++;
        else
          finished = false;
      }
    }
    else
    {
      if (!depth)
      {
        if (finished)
          break;
        status(Priority_General,2,1,"Examining word-acceptor to depth %d\n",++limit);
        we.lhs_word().set_length(0);
        finished = true;
      }
      else
        depth--;
    }
  }
  delete [] state;
  if (finished)
    container.progress(1,"Exploration of word-acceptor completed\n");
  else
    container.progress(1,"Exploration of word-acceptor timed out\n");
  update_machine(UM_CHECK_POOL);
  update_machine(UM_RECHECK_PARTIAL+UM_NO_REPEAT_PARTIAL);
  State_Count d;
  while ((d = dt->take_changes())!=0)
  {
    stats.difference_delta += d;
    update_machine(UM_RECHECK_PARTIAL+UM_NO_REPEAT_PARTIAL);
    retcode = true;
  }
  return retcode;
}

/**/

void Rewriter_Machine::explore_dm(bool repeat)
{
  /* This method uses the current difference tracker to reduce the
     words corresponding to nodes in the Rewriter_Machine. We do this
     because spurious word-differences may have crept in because of
     incomplete reductions.
     The "repeat" parameter will be true when we want to build an
     automatic structure right now, and false if we just want to
     remove as many spurious word-differences as we can.
  */

  if (!dt || !pd.is_short && !repeat)
    return;
  stats.last_explore = stats.status_count;

  for (;;)
  {
    dt->scan_all_differences();
    Total_Length vital_limit = dt->vital_eqn_limit();
    if (stats.visible_limit < vital_limit && vital_limit <= 12)
      stats.visible_limit = vital_limit;
    update_machine();

    dt->compute_transitions(0);
    State_ID nr_states = dt->state_count();
    Strong_Diff_Reduce sdr(this);
    Node_Count added = 0;
    Ordinal_Word lhs_word(alphabet(),1);
    Ordinal_Word rhs_word(alphabet(),1);
    Ordinal nr_g_generators = dt->base_alphabet.letter_count();
    Equation_Word ew2(nm);
    Equation_Word ew3(nm);
    Working_Equation we(nm,Derivation(BDT_Diff_Reduction));
    /* We should have an option to disable the checks done in this for
       loop, since they can take a long time and not achieve anything useful
    */
    for (State_ID state = 2; state < nr_states && !maf.aborting;state++)
    {
      Node_Reference s = dt->nm_state(state);
      if (s.is_null())
        continue;
      s->read_word(&lhs_word,nm);
      status(Priority_General,
             2,1,"Reducing words near word-differences ("
                 FMT_ID " of " FMT_ID ") (" FMT_NC " reductions found)\n",
             state,nr_states,added);
      Word_Length l = lhs_word.length();
      rhs_word.set_length(1);
      rhs_word.set_code(0,0);
      lhs_word = rhs_word + lhs_word;
      for (Ordinal g1 = 0;g1 < nr_g_generators;g1++)
      {
        if (!redundant(inverse(g1)))
        {
          lhs_word.set_length(l+1);
          lhs_word.set_code(0,g1);
          ew2 = lhs_word;
          ew2.reduce();
          for (Ordinal g2 = 0;g2 < nr_g_generators;g2++)
          {
            if (!redundant(g2))
            {
              ew3 = ew2;
              ew3.append(g2);
              ew3.reduce();
              s = ew3.state_word(true);
              if (s.is_null())
              {
                if (repeat && sdr.reduce(&rhs_word,ew3))
                {
                  we.lhs_word() = ew3;
                  we.rhs_word() = rhs_word;
                  add_equation(&we,AE_KEEP_LHS|AE_INSERT|AE_NO_REPEAT);
                  added++;
                }
              }
              else
                s->node(nm).set_flags(NF_NEAR_DIFFERENCE);
            }
          }
        }
      }
      if (pd.is_normal_coset_system && dt->is_initial(state) && !stats.coset_complete)
      {
        /* The initial states of the tracker generate the subgroup, so check that
           all conjugates of the generators are in the subgroup */
        s = dt->nm_state(state);
        if (s.is_null())
          continue;
        if (dt->nr_initial_states() > alphabet().capacity())
          MAF_INTERNAL_ERROR(container,
                             ("Normal subgroup has too many generators, and"
                              " may not be finitely generated\n"));
        s->read_word(&lhs_word,nm);
        rhs_word.set_length(0);
        rhs_word.append(pd.coset_symbol);
        rhs_word.append(0);
        lhs_word = rhs_word + lhs_word;
        for (Ordinal g = 0;g < nr_g_generators;g++)
        {
          rhs_word.set_code(1,g);
          lhs_word.set_code(1,g);
          Working_Equation we(nm,lhs_word,rhs_word,Derivation(BDT_Normal_Closure,s,g));
          if (add_equation(&we,AE_KEEP_LHS|AE_INSERT|AE_NO_REPEAT))
            added++;
        }
      }
    }
    update_machine();
    if (dt->state_count() != nr_states)
      continue;

    Word_Length depth = 0;
    Node_Count done = 0;
    Equation_Word ew(nm);
    ew.set_length(0);
    ew.reduce();
    ew2.set_length(0);
    ew2.reduce();

    added = 0;
    bool do_all = !repeat;
    bool restart = true;
    Ordinal child_start,child_end;
    while (!maf.aborting)
    {
      /* In this block we look for nodes that can be reduced using the
         the Difference_Tracker. If repeat is false, then initially
         we only look for reductions that definitely affect the word
         differences:

         1) LHS prefixes and RHS of equations that "own" a word-difference. If
          we find a reduction in one of these nodes, we will switch to doing a
          full pass of the tree.

          2) Word-differences, half-differences and nodes near word-differences.
          Finding a reduction here is not quite so bad: we may be able to
          relabel the difference, but will do a full pass if we find the word
          difference machine has changed.
      */

      if (restart)
      {
        done = added = 0;
        depth = 0;
        ew.set_length(0);
        restart = false;
      }

      Ordinal v = ew.values[depth]+1;
      nm.valid_children(&child_start,&child_end,depth ? ew.values[depth-1] : IdWord);
      if (v < child_start)
        v = child_start;
      if (v < child_end)
      {
        ew.set_length(depth);
        ew.append(v);
        /* See comment above in explore_acceptor() */
        Node_Reference ns = ew.state_lhs(false);
        if (ns.is_null())
        {
          if (ew.reduce(AE_KEEP_LHS))
          {
            /* In this case we have managed to delete the part of
               the tree we were working on, so we must have created an
               equation and need to start again */
            restart = true;
            continue;
          }
          ns = ew.state_lhs(false);
        }

        if (!ns.is_null() && ns->length(nm) == depth+1)
        {
          status(Priority_General,2,1,
                 "Reducing nodes with difference machine (" FMT_NC " of " FMT_NC ") "
                 "(" FMT_NC " equations added)\n",done++,
               nm.stats.nc[0]+nm.stats.nc[1]+nm.stats.nc[2]+nm.stats.nc[3],added);
          int do_it = do_all;
          if (!ns->is_final())
          {
            if (!do_it)
            {
              do_it = ns->flagged(NF_IS_DIFFERENCE+NF_IS_HALF_DIFFERENCE+NF_NEAR_DIFFERENCE);
              for (State e = ns->first_lhs(nm);e; e = e->next_lhs(nm))
                if (e->flagged(EQ_HAS_DIFFERENCES))
                {
                  do_it = 2;
                  break;
                }
            }
            depth++;
          }
          else
          {
            do_it = ns->flagged(EQ_HAS_DIFFERENCES) ? 2 : 0;
            if (ns->fast_is_primary() && (do_it || !repeat))
            {
              ns->attach(nm);
              for (Word_Length i = 1; i < depth;i++)
              {
                lhs_word = Subword(ew,i);
                rhs_word = lhs_word;

                if (sdr.reduce(&rhs_word,rhs_word,DR_ALTERNATE))
                {
                  we.lhs_word() = lhs_word;
                  we.rhs_word() = rhs_word;
                  if (add_equation(&we,AE_KEEP_LHS|AE_INSERT|AE_NO_REPEAT)==1)
                  {
                    update_machine();
                    added++;
                  }
                }
                else
                  break;
              }
              if (ns->flagged(NF_REMOVED))
                do_it = false;
              ns->detach(nm,ns);
            }
          }

          if (do_it)
          {
            Ordinal_Word rhs_word(ew);
            if (ns->is_final())
              rhs_word.set_length(depth);
            if (sdr.reduce(&rhs_word,rhs_word))
            {
              we.lhs_word() = ew;
              if (ns->is_final())
                we.lhs_word().set_length(depth);
              we.rhs_word() = rhs_word;
              if (!ns->is_final())
                depth--;
              if (add_equation(&we,AE_KEEP_LHS|AE_INSERT|AE_NO_REPEAT)==1)
              {
                update_machine();
                added++;
                if (do_it == 2)
                {
                  restart = do_all = true;
                  continue;
                }
              }
            }
          }
        }
      }
      else
      {
        if (!depth)
          break;
        else
          depth--;
      }
    }
    if (!repeat || maf.aborting || !do_all && dt->state_count()==nr_states)
      break;
  }
  stats.explore_time = stats.status_count - stats.last_explore;
  stats.last_explore = stats.status_count;
}

/**/

void Rewriter_Machine::schedule_right_conjugation(Equation_Handle e)
{
  if (e->is_valid())
  {
    Ordinal ig = inverse(e->last_letter());
    Node_Reference rhs_node = e->raw_reduced_node(nm);
    Node_Reference prefix = e->parent(nm);
    if (ig != INVALID_SYMBOL)
    {
      Node_Reference conjugate = rhs_node->transition(nm,ig);
      if (conjugate->reduced_node(nm,conjugate) != prefix &&
          rhs_node->length(nm) < MAX_WORD)
      {
        Derivation d(BDT_Right_Conjugate,e);
        if (proof_trace)
          d.set_trace(proof_trace->node_equation(nm,e),0);
        jm.equate_nodes(rhs_node,prefix,d,ig,true);
      }
    }
  }
}

/**/

bool Rewriter_Machine::conjugate(Equation_Handle e)
{
  /* This function changes the equation ax=b to the  equation x=Ab.
     If the full parameter is true, which it will be unless very many
     equations need conjugating at the same time, overlaps of which the LHS
     is the prefix or trailing subword are processed.
  */
  if (maf.aborting)
    return false;
  Node_Reference rhs_node = e->raw_reduced_node(nm);
  Word_Length l = e->lhs_length(nm);
  if (maf.options.conjugation >= 1 &&
      !(pd.is_coset_system && maf.options.no_h && e->last_letter() > pd.coset_symbol))
  {
    Word_Length r = rhs_node->irreducible_length();
    Ordinal lvalue = e->first_letter(nm);
    Ordinal ilvalue = inverse(lvalue);
    if (ilvalue != INVALID_SYMBOL) // && l+r <= stats.discard_limit)
    {
      Ordinal_Word lhs_word(alphabet(),l);
      Ordinal_Word rhs_word(alphabet(),r+1);
      Ordinal * rvalues = rhs_word.buffer();
      e->read_values(lhs_word,nm,l);
      rhs_node->read_values(rvalues+1,nm,r);
      rvalues[0] = ilvalue;
      lhs_word = Subword(lhs_word,1,WHOLE_WORD);
      int primary = e->fast_is_primary();
      Working_Equation we(nm, primary ? rhs_word : lhs_word,
                          primary ? lhs_word : rhs_word,
                          Derivation(BDT_Left_Conjugate,e,ilvalue));
      unsigned ae_flags = AE_DISCARDABLE|AE_INSERT_DESIRABLE;
      if (stats.want_secondaries || !primary)
        ae_flags |= AE_KEEP_LHS;
      add_equation(&we,ae_flags);
    }

    if (maf.options.conjugation == 4 && lvalue < pd.coset_symbol)
    {
      Ordinal_Word lhs_word(alphabet(),l+1);
      Ordinal_Word rhs_word(alphabet(),r+1);
      Ordinal * lvalues = lhs_word.buffer();
      Ordinal * rvalues = rhs_word.buffer();
      e->read_values(lvalues+1,nm,l);
      rhs_node->read_values(rvalues+1,nm,r);
      lvalues[0] = rvalues[0] = pd.coset_symbol;
      Working_Equation we(nm, lhs_word,rhs_word,
                          Derivation(BDT_Left_Conjugate,e,pd.coset_symbol));
      unsigned ae_flags = AE_DISCARDABLE|AE_INSERT_DESIRABLE;
      add_equation(&we,ae_flags);
    }

    if (maf.options.conjugation>=2 &&
        l+r+2 <= stats.visible_limit && !nm.is_confluent() && e->is_valid() &&
        maf.options.conjugation != 4)
    {
      Simple_Equation se1(nm,e);

      if (e->fast_is_primary() && l != 1)
      {
        Node_Reference trailing_subword = e->overlap_suffix(nm);
        Ordinal g = e->last_letter();
        Ordinal ig = inverse(g);
        trailing_subword->attach(nm);
        Ordinal child_start,child_end;
        pd.right_multipliers(&child_start,&child_end,g);
        for (g = child_start; g < child_end; g++)
        {
          if (g != ig) /* we already formed the right conjugate */
          {
            Node_Reference e2 = trailing_subword->transition(nm,g);
            if (e2->is_equation() &&
                (maf.options.conjugation==3 ?
                 e2->lhs_length(nm) > 1 : e2->lhs_length(nm) == 2))
            {
              Word_Length offset = l+1-e2->lhs_length(nm);
              Word_Length nl = offset + e2->reduced_length(nm);
              Word_Length nr = r + 1;
              if (nl+nr <= l+r || nl + nr <= e2->total_length(nm) ||
                  nl+nr <= stats.visible_limit)
              {
                Simple_Equation se2(nm,e2);
                e2->attach(nm);
                consider(e,e2,offset,se1,se2,AE_DISCARDABLE);
                e2->detach(nm,e2);
                if (!e->is_valid() || !e->fast_is_primary())
                  break;
              }
            }
          }
        }
        trailing_subword->detach(nm,trailing_subword);
      }

      if (e->is_valid() && lvalue != pd.coset_symbol)
      {
        Ordinal left_start,left_end;
        pd.left_multipliers(&left_start,&left_end,lvalue);
        Word_Length limit = maf.options.conjugation==3 ? l-1 : 1;
        for (Ordinal g = left_start; g < left_end; g++)
        {
          State s = nm.start()->transition(nm,g);
          if (g != ilvalue && !s->is_final())
          {
            Word_Length i;
            for (i = 0; i < limit ;i++)
            {
              Node_Reference e2 = s->transition(nm,se1.lhs_word.value(i));
              if (e2->is_equation())
              {
                Word_Length nl = e2->reduced_length(nm) + l - i;
                Word_Length nr = 1 + r;
                if (nl+nr <= l+r || nl + nr <= e2->total_length(nm) ||
                    nl+nr <= stats.visible_limit)
                {
                  Ordinal_Word lhs_word(alphabet(),l+1);
                  Ordinal_Word rhs_word(alphabet(),r+1);
                  rhs_word.set_length(1);
                  rhs_word.set_code(0,g);
                  lhs_word = rhs_word + se1.lhs_word;
                  rhs_word += se1.rhs_word;
                  Working_Equation we(nm,rhs_word,lhs_word,
                                      Derivation(BDT_Left_Conjugate,e,g));
                  add_equation(&we,AE_DISCARDABLE);
                }
                if (!e->is_valid())
                  return false;
                break;
              }
            }
          }
        }
      }
    }
  }
  if (maf.options.check_inverses)
    rhs_node->node(nm).set_inverse(nm,rhs_node);
  return true;
}

/**/

void Rewriter_Machine::improve_rhs(Equation_Handle e)
{
  /* Improve the RHS of an equation that was on the list of equations with
     a bad RHS. It may be that the new RHS has a common prefix or suffix
     with the LHS, or that balancing shows that a prefix of the LHS is
     reducible. In that case the equation is now bad and should not be in
     the tree at all. This is the reason we mess around with the
     AE_MUST_SIMPLIFY and EQ_REDUNDANT flags below.

     We also have to deal with the unpleasant situation where the new
     RHS is too long, and so the equation cannot be repaired.

     Even if an equation being repaired keeps the same LHS it is
     a new equation from the point of view of later processing.
  */

  Working_Equation we(nm,e,
                      maf.options.log_flags & LOG_DERIVATIONS ? e->raw_reduced_node(nm) :
                      e->fast_reduced_node(nm),Derivation(BDT_Known,e));
  unsigned ae_flags = AE_COPY_LHS|AE_KEEP_LHS|AE_NO_UNBALANCE|AE_INSERT|AE_NO_REPEAT|AE_KEEP_FAILED;
  /* See if there might be a prefix change link pointing to this equation.
     If there could be then we will use check_transitions() to ensure it
     gets removed if it is redundant. If not then we can get rid of it
     cheaply by replacing it with its suffix. */
  if (e->reference_count > 2 && e->flagged(NF_TARGET))
    ae_flags |= AE_MUST_SIMPLIFY;
  for (;;)
  {
    /* We set EQ_REDUNDANT here, and check if it is still set after we
       have added the equation. If it is then either the equation is
       redundant but we did not set the AE_MUST_SIMPLIFY flag,
       or it is not redundant at all, but we could not repair the equation
    */
    e->node(nm).set_flags(EQ_REDUNDANT);
    add_equation(&we,ae_flags);
    if (e->flagged(NF_REMOVED) || !e->flagged(EQ_REDUNDANT))
      break; /* We either repaired the equation or this part of the tree has now gone altogether */

    /* If we get here then the equation is still in the tree, but
       ought not to be */
    we.lhs_word().read_word(e,false);
    Node_Reference ns = nm.get_state(((const Equation_Word &) we.lhs_word()).buffer(),
                                      e->length(nm),1);
    if (ns->is_final())
    {
      /* replace the equation by a prefix change */
      Node_Reference prefix = e->parent(nm);
      Ordinal rvalue = e->last_letter();
      prefix->node(nm).link_replace(nm,prefix,rvalue,ns,e);
      jm.schedule_partial_reduction_check(prefix,rvalue);
      if (ae_flags & AE_MUST_SIMPLIFY)
      {
        we.spare->read_word(ns,true);
        nm.check_transitions(*we.spare,we.spare->word_length,false);
      }
      break;
    }
    else
    {
      /* The equation cannot be repaired! */
      e->node(nm).clear_flags(EQ_REDUNDANT);
      stats.complete = false;
      break;
    }
  }
#if 0
  if (e->is_valid())
  {
    we.derivation.reset(Derivation(BDT_Known,e));
    we.changed = false;
    if (we.normalise(0) && we.changed)
      add_equation(&we,AE_NO_REPEAT+AE_DISCARDABLE);
  }
#endif
}

/**/

void Rewriter_Machine::improve_dubious(Equation_Handle e)
{
  /* Check the RHS of an equation that was on the list of equations with
     a dubious RHS.
  */

  /* If the equation is valid put it back on the list it was on before -
     consider will remove it if necesary */
  if (!e->raw_reduced_node(nm)->is_final())
  {
    jm.cancel_jobs(e);
    switch (Partial_Switch(e->status()))
    {
      case NS_Dubious_Secondary_Equation:
        jm.schedule_jobs(e,NS_Weak_Secondary_Equation);
        break;
      case NS_Dubious_Unconjugated_Equation:
        jm.schedule_jobs(e,NS_Unconjugated_Equation);
        break;
      case NS_Dubious_Undifferenced_Equation:
        jm.schedule_jobs(e,choose_status(e,NS_Undifferenced_Equation));
        break;
    }
  }

  Node_Reference primary = e->fast_primary(nm);
  Simple_Equation se(nm,e);
  Simple_Equation se2(nm,primary);
  Word_Length offset = se.lhs_length-se2.lhs_length;
  Word_Length r = offset+se2.rhs_length;
  Ordinal_Word rhs_word(alphabet(),r);
  word_copy(rhs_word,se.lhs_word,offset);
  word_copy(rhs_word.data(offset),se2.rhs_word,se2.rhs_length);
  Working_Equation we(nm,se.rhs_word,rhs_word,
                      Derivation(BDT_Conjunction,e,primary,offset));
  bool better = false;
  bool conflict = false;
  bool resolved = false;

  if (we.normalise(AE_REDUCE|AE_SAFE_REDUCE|AE_NO_BALANCE|AE_PRE_CANCEL))
  {
    better = we.swapped;
    conflict = true;
    /* AE_INSERT here is not safe */
    resolved = add_equation(&we,AE_NO_UNBALANCE|AE_NO_REPEAT)==1;
    if (resolved && se.rhs_length != 0 && maf.options.weed_secondary)
      better = false;
  }
  else
    better = e->reduced_length(nm)==0;// This is a surprisingly common case
                                    // when recursive type orders are used

  if (maf.options.secondary==0)
  {
    better = false;
    conflict = true;
  }

  if (!e->flagged(NF_REMOVED))
  {
    if (stats.want_secondaries || better)
    {
      if (e->raw_reduced_node(nm)->is_final())
        improve_rhs(e);
      else if (better && !resolved &&
               e->status() == NS_Weak_Secondary_Equation &&
               maf.options.consider_secondary)
      {
        jm.cancel_jobs(e);
        jm.schedule_jobs(e,NS_Adopted_Equation);
      }
    }
    else if (conflict || !e->is_valid() || maf.options.weed_secondary)
    {
      if (!pd.is_short && !resolved)
      {
        bool failed = false;
        Equation_Word ew(nm,rhs_word);
        ew.reduce(0,&failed);
        if (!failed && ew.word_length < height() && ew.word_length*2 > stats.visible_limit)
          ew.realise_state();
      }

      // The primary of e may have changed but the old value must still be
      // a valid equation, so it makes sense to use this to get rid of
      // e. Selecting either the longest remaining secondary or the
      // new primary does not seem to improve matters at all.
if (primary->flagged(NF_REMOVED))
 * (char *) 0 = 0;
      if (e->reference_count > 2 && e->flagged(NF_TARGET))
      {
        Equation_Word ew(nm,primary,-1);
        e->node(nm).set_flags(EQ_REDUNDANT);
        nm.check_transitions(ew,ew.word_length,false);
      }
      else
      {
        Node_Reference prefix = e->parent(nm);
        e->node(nm).set_flags(EQ_REDUNDANT);
        Ordinal rvalue = e->last_letter();
        prefix->node(nm).link_replace(nm,prefix,rvalue,primary,e);
        jm.schedule_partial_reduction_check(prefix,rvalue);
      }
    }
  }
}

/**/

Node_Status Rewriter_Machine::learn_differences(Equation_Handle e)
{
  /* For a reducible node e, recreate the equation and teach it to
     Difference_Tracker */
  if (dt)
  {
    Working_Equation we(nm,e,e->fast_reduced_node(nm),Derivation(BDT_Known,e));
    Total_Length total = we.total_length();
    if (!pd.is_coset_system || !e->parent(nm)->flagged(NF_INVISIBLE))
    {
      if (total > dt->status().max_difference_eqn &&
          total > stats.visible_limit && total > 2 &&
          !dt->is_accepted_equation(we) && !pd.inversion_difficult)
      {
        /* Here we have a long equation which has survived conjugation and is
           about to impact the word-differences. Before it does so let us make
           sure it is not spurious. This also helps to stabilise the word
           differences (e.g. in presentations like "trefoil" or "cox33334c"
           where otherwise spurious word-differences recur for a long time */
        Diff_Reduce dr(dt);
        Ordinal_Word lhs_word(we.lhs_word());
        Ordinal_Word rhs_word(we.rhs_word());
        if (dr.reduce(&lhs_word,lhs_word,DR_PREFIX_ONLY))
        {
          Working_Equation we2(nm,we.lhs_word(),lhs_word,
                              Derivation(BDT_Diff_Reduction));
          add_equation(&we2,AE_KEEP_LHS|AE_INSERT|AE_NO_REPEAT);
        }

        if (dr.reduce(&rhs_word,rhs_word))
        {
          Working_Equation we2(nm,we.rhs_word(),rhs_word,
                              Derivation(BDT_Diff_Reduction));
          add_equation(&we2,AE_KEEP_LHS|AE_INSERT|AE_NO_REPEAT);
        }
        if (!e->is_valid())
          return NS_Weak_Secondary_Equation;
      }
      dt->set_limit(stats.visible_limit);
      if (dt->learn_equation(we,e,PD_CHECK))
        stats.last_difference = stats.status_count;
    }
  }
  return choose_status(e,maf.options.strategy & MSF_CONJUGATE_FIRST ?
                       NS_Adopted_Equation : NS_Unconjugated_Equation);
}

/**/

unsigned Rewriter_Machine::learn_differences(Working_Equation & we)
{
  /* Teach an equation to Difference_Tracker */
  return dt->learn_equation(we,Node_Reference(0,0),PD_CHECK);
}

/**/

Node_Status Rewriter_Machine::choose_status(Equation_Handle e,
                                            Node_Status recommendation)
{
  if (e->is_valid())
  {
    Total_Length tl = e->total_length(nm);
    if (recommendation == NS_Undifferenced_Equation)
    {
      if (stats.want_differences || !stats.started)
      {
        if (tl > stats.visible_limit*2 && pd.inversion_difficult &&
            !nm.is_confluent() && !pd.is_coset_finite)
          return NS_Oversized_Equation;
        return recommendation;
      }
      if (!(maf.options.strategy & MSF_CONJUGATE_FIRST))
        return NS_Unconjugated_Equation;
    }
    if (recommendation == NS_Unconjugated_Equation)
      return NS_Unconjugated_Equation;

    if (e->fast_is_primary())
    {
      Word_Length l = e->lhs_length(nm);
      if ((maf.options.no_pool || maf.options.expand_all) &&
          tl > stats.max_expanded || l > stats.max_expanded_lhs)
        return NS_Adopted_Equation;

      if (l >= stats.max_primary_lhs)
        stats.max_primary_lhs = l;
      if (tl <= stats.auto_expand_limit ||
          pd.is_coset_system && maf.options.detect_finite_index==2 &&
          e->first_letter(nm) == pd.coset_symbol)
      {
        if (tl > stats.max_expanded)
          stats.max_expanded = tl;
        if (l > stats.max_expanded_lhs)
          stats.max_expanded_lhs = l;
        stats.approved_complete = false;
        return NS_Expanded_Equation;
      }

      if (!stats.pool_full)
      {
        if (!pd.g_level)
        {
          if (tl <= stats.max_expanded &&
              (e->length(nm) <= stats.lhs_limit &&
               2*e->reduced_length(nm) <= stats.auto_expand_limit))
            return NS_Expanded_Equation;
        }
        if (!pd.inversion_difficult)
        {
          if (stats.favour_differences && stats.vital_attained && e->flagged(EQ_INTERESTING) &&
              (tl <= stats.max_expanded || tl <= stats.visible_limit))
          {
            if (l > stats.max_expanded_lhs)
              stats.max_expanded_lhs = l;
            if (e->flagged(EQ_HAS_DIFFERENCES+EQ_HAS_PRIMARY_DIFFERENCES))
              return NS_Expanded_Equation;
            if (tl <= stats.max_expanded && tl <= stats.visible_limit)
              return NS_Expanded_Equation;
          }
          if (stats.favour_differences && e->flagged(EQ_HAS_DIFFERENCES+EQ_HAS_PRIMARY_DIFFERENCES))
            if (l > stats.max_expanded_lhs)
              stats.max_expanded_lhs = l;
        }
      }
    }
    else
      return stats.want_overlaps ? NS_Dubious_Secondary_Equation :
                                   NS_Weak_Secondary_Equation;
  }
  return NS_Adopted_Equation;
}

/**/

bool Rewriter_Machine::relabel_difference(Node_Handle node,Node_Flags old_flags,Node_Count to_do)
{
  bool retcode = false;
  if (to_do > 100 && to_do > height()*jm.status().nr_visible)
  {
    /* If there are more difference corrections pending than there can
       possibly be word-differences we should delete the word-difference machine
       because it is useless.
       This will usually happen after a collapse and the required amalgamation
       of states can take a very long time
    */
    delete dt;
    dt = 0;
  }
  if (dt)
  {
    if (dt->relabel_difference(node,old_flags))
      stats.last_difference = stats.status_count;
    retcode = true;
  }
  if (!to_do && !dt)
  {
    /* we have no real idea how many many word-differences exist now
       so revert to the initial default of 1024 */
    dt = Difference_Tracker::create(nm,1024);
  }
  return retcode;
}

/**/

Word_Length Rewriter_Machine::height() const
{
  return nm.height();
}

/**/

Node_Count Rewriter_Machine::node_count(bool include_equations) const
{
  if (include_equations)
    return nm.stats.nc[language_L0] + nm.stats.nc[language_L1] +
           nm.stats.nc[language_L2] + nm.stats.nc[language_L3];

  return nm.stats.nc[language_L0];
}

/**/

Node_Count Rewriter_Machine::primary_count() const
{
  return nm.stats.nc[language_L1];
}

/**/

void Rewriter_Machine::prune()
{
  if (stats.favour_differences && pd.is_short && !maf.aborting && !maf.options.no_prune)
  {
    Word_Length prune_length = stats.visible_limit/2+1;

    if (prune_length < dt->status().max_difference_length)
    {
      check_differences(0);
      if (dt->status().max_difference_length > prune_length)
        prune_length = dt->status().max_difference_length;
    }

    update_machine();
    Node_Iterator ni(nm);
    Node_Reference s;
    Node_List prune_list;
    Node_Count to_do = 0;
    bool inside = true;
    while (ni.next(&s,inside))
    {
      inside = ni.word().length() < prune_length;
      if (!inside)
        if (!s->is_final() || s->total_length(nm) > stats.visible_limit)
          prune_list.add_state(s);
      to_do = prune_list.length();
      container_status(2,1,"Building list of long nodes. " FMT_NC " found so far\n",to_do);
    }

    Node_Count done = 0;

    while (prune_list.use(&s,nm))
    {
      Equation_Word ew(nm,s);
      container_status(2,1,"Pruning " FMT_NC " of " FMT_NC " long nodes\n",done++,to_do);
      Node_Reference ns = nm.get_state(ew.buffer(),s->length(nm),1);
      Node_Reference parent = s->parent(nm);
      parent->node(nm).link_replace(nm,parent,s->last_letter(),ns,s);
    }
  }
}

/**/

bool Rewriter_Machine::need_partial_reduction_checks() const
{
  if (!dt)
    return false;
  if (maf.options.assume_confluent && stats.complete && stats.want_differences)
    return true;
  const Difference_Tracker::Status & dstats = dt->status();
  return dstats.nr_tertiary_differences != 0;
}

/**/

Total_Length Rewriter_Machine::interest_limit() const
{
  if (!dt)
    return 0;
  return dt->interest_limit();
}

/**/

bool Rewriter_Machine::dm_valid(unsigned gwd_flags)
{
  return dt ? !dt->dm_changed(gwd_flags,false) : true;
}

/**/

bool Rewriter_Machine::is_dm_changing(bool include_history,unsigned gwd_flags,bool force_tick)
{
  if (!dt)
    return true;
  State_Count d = dt->take_changes();
  bool changed = (d != 0);
  if (stats.normal_timeout < maf.options.timeout/2 && d > 100)
    stats.normal_timeout++;
  stats.difference_delta += d;
  if (stats.status_count - stats.ticked >= stats.timeout || force_tick)
  {
    stats.difference_delta /= 2;
    if (stats.status_count - stats.last_difference >= maf.options.timeout &&
        stats.status_count - stats.expand_time >= maf.options.timeout)
      stats.difference_delta = 0;
    stats.timeout = stats.normal_timeout;
    stats.ticked = stats.status_count;
    stats.just_ticked = true;
    if (stats.normal_timeout < maf.options.timeout/2 && dt->take_max_changed())
      stats.normal_timeout++;
  }
  if (gwd_flags)
    changed = dt->take_changes(gwd_flags)!=0;
  return changed || include_history && stats.difference_delta != 0;
}

/**/

State_ID Rewriter_Machine::get_derivation_id(Node_ID eid, bool insert)
{
  if (!derivation_db)
  {
    derivation_db = new Hash(1024*1024,sizeof(Node_ID));
    const Equation_Node * zero = 0; // Insert dummy equation so we start from 1
    derivation_db->find_entry(&zero,sizeof(Node_ID));
  }
  State_ID answer = derivation_db->find_entry(&eid,sizeof(Node_ID),insert);
  return answer > 0 ? answer : 0;
}

/**/

Strong_Diff_Reduce::Strong_Diff_Reduce(Rewriter_Machine * rm_) :
  rm(rm_),
  dt(rm_ ? rm_->dt : 0)
{
  dr0 = dt ? new Diff_Reduce(dt) : 0;
  dr1 = dt ? new Diff_Reduce(dt) : 0;
  dr2 = 0;
}

Strong_Diff_Reduce::~Strong_Diff_Reduce()
{
  if (dr0)
    delete dr0;
  if (dr1)
    delete dr1;
  if (dr2)
    delete dr2;
}

/**/

unsigned Strong_Diff_Reduce::reduce(Word * rword,const Word & word,
                                    unsigned flags,const FSA *)
{
  if (rm->stats.complete && !(flags & WR_ONCE))
    return rm->reduce(rword,word,flags);
  if (dt->dm_changed(0,true))
  {
    dr0->invalidate();
    dr1->invalidate();
    if (dr2)
      dr2->invalidate();
    dt->clear_changes();
  }
  int retcode = 0;
  if (flags & DR_ALTERNATE)
  {
    if (!dr2)
      dr2 = new Diff_Reduce(dt);
    if (dr2->reduce(rword,word,flags|DR_ONCE|DR_NO_G_PREFIX))
      retcode = 1;
  }
  else if (dr0->reduce(rword,word,flags|DR_ONCE|DR_NO_G_PREFIX))
    retcode = 1;
  if (retcode)
    retcode += dr1->reduce(rword,*rword,flags|DR_NO_G_PREFIX);
  return retcode;
}

/**/

void Rewriter_Machine::import_difference_machine(const FSA_Simple & dm)
{
  /* This function allows us to import the states and transitions of
     an externally built word-difference machine into a Rewriter_Machine
     instance. This is for times when MAF is not trying to discover an
     automatic structure, but work with one that has previously been built
  */
  for (;;)
  {
    clear_differences(0);
    if (maf.options.no_equations)
      maf.options.no_half_differences = true;
    if (!dt->import_difference_machine(dm))
      break;
    update_machine();
    explore_dm(false);
    if (dt->is_valid(0) && dt->is_inverse_complete())
      break;
    check_differences(CD_STABILISE);
  }
}

/**/


Total_Length Rewriter_Machine::overlap_size(State e1,
                                            State e2,
                                            Word_Length offset) const
{
  return e1->reduced_length(nm) + e2->length(nm) - e1->length(nm) + offset +
         offset + e2->reduced_length(nm);
}

bool Rewriter_Machine::overlap_interesting(State e1,
                                           State e2,
                                           Word_Length offset)
{
  if (!maf.options.special_overlaps)
    return false;
  if (maf.options.no_h && e2->first_letter(nm) > pd.coset_symbol)
    return false;
  if (maf.options.no_pool)
  {
    if (e1->flagged(EQ_EXPANDED) || !e1->expanded_timestamp())
      return false;
    if (maf.options.strategy & MSF_LEFT_EXPAND || stats.starting)
      return false;
    /* there used to be a return true here, but that makes no_pool very
       unsafe unless special_overlaps are disabled completely */
  }
  Total_Length tl = overlap_size(e1,e2,offset);
  if (tl <= stats.special_limit)
    return true;
  if (maf.options.special_overlaps <= 2)
    return false;
  return tl <= e1->total_length(nm) || tl <= e2->total_length(nm) ||
         maf.options.special_overlaps >= 4;
}
