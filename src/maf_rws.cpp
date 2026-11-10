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


// $Log: maf_rws.cpp $
// Revision 1.15  2010/06/10 13:57:40Z  Alun
// All tabs removed again
// Revision 1.14  2010/04/10 19:11:45Z  Alun
// Changes required by new style Node_Manager interface
// Progress messages improved
// Revision 1.13  2009/09/16 07:39:55Z  Alun
// Additional source code changes needed to get clean compilation on latest GNU
//
// Revision 1.12  2009/09/13 19:05:52Z  Alun
// printf() style functions cleaned up considerably
// code now compiles cleanly in g++ both 32-bit and 64-bit
// all tabs removed from source files
// Revision 1.11  2008/12/28 21:22:58Z  Alun
// Changed management of NF_NEEDED flag because it now has same value as NF_NEW
// Revision 1.13  2008/11/05 23:35:05Z  Alun
// Last version had a typo that prevented .kbprog files from being created properly
// Revision 1.12  2008/11/04 23:32:27Z  Alun
// Added facility to use RWS constructor as way of building word-acceptor
// Revision 1.11  2008/11/02 18:57:14Z  Alun
// word reduction improved for recursive orderings as per maf_we.cpp code
// Revision 1.10  2008/10/01 06:32:55Z  Alun
// Change to Hash broke code as I forgot to ask for eqn to be managed
// Revision 1.9  2008/09/29 19:58:04Z  Alun
// Switch to using Hash+Collection_Manager as replacement for Indexer.
// Currently this is about 10% slower, but this is a more flexible
// architecture.
// Revision 1.4  2007/12/20 23:25:42Z  Alun
//

/* This module contains code to create a KBMAG style RWS from a
   Rewriter_Machine, and to save and load such objects to/from a KBMAG
   format file.
*/

#include "maf.h"
#include "mafword.h"
#include "keyedfsa.h"
#include "maf_rws.h"
#include "maf_rm.h"
#include "mafnode.h"
#include "maf_nm.h"
#include "container.h"
#include "maf_trace.h"

Rewriting_System::~Rewriting_System()
{
  delete &eqn;
  delete [] inverses;
}

/**/

Rewriting_System::Rewriting_System(Rewriter_Machine * rm,unsigned flags,
                                   FSA_Simple **taker) :
  Delegated_FSA(rm->container,rm->alphabet(),true),
  eqn(* new Hash(rm->primary_count()+1,sizeof(Node_ID))),
  source(flags & RWSC_MINIMAL ? rm : 0)
{
  /* this function builds a KBMAG style rewriting system from a MAF
     Rewriter_Machine() */
  Ordinal nr_generators = rm->generator_count();
  const Presentation_Data & pd = rm->pd;
  bool want_nice_coset_system = pd.is_coset_system;
  Node_Manager & nm = *rm;
  if (flags & RWSC_NO_H)
    nr_generators = pd.coset_symbol+1;

  /* Set up miscellaneous information needed for correct output */
  inverses = new Ordinal[nr_generators];
  for (Ordinal i = 0; i < nr_generators;i++)
    inverses[i] = rm->inverse(i);
  confluent = nm.is_confluent();

  Ordinal * transition_order = new Ordinal[nr_generators];
  Ordinal i;
  /* An alternative state ordering might be nice for other non shortlex
     systems as well. One way to get it would be to put the generators
     into an Ordered_Node_Tree and pluck them out. In fact one could
     even push all the states that will be needed into the tree */

  if (want_nice_coset_system)
  {
    int j = 0;
    transition_order[j++] = pd.coset_symbol;
    for (i = 0;i < pd.coset_symbol; i++)
      transition_order[j++] = i;
    for (i = pd.coset_symbol+1; i < nr_generators;i++)
      transition_order[j++] = i;
  }
  else
    for (i = 0; i < nr_generators;i++)
      transition_order[i] = i;

  eqn.manage(packed_lhses);
  eqn.manage(packed_rhses);
  /* Now create the RWS fsa and equation set */
  State_ID * transition = new State_ID[nr_generators];
  Keyed_FSA factory(container,rm->alphabet(),nr_generators,
                             rm->node_count()+1,sizeof(Node_ID),false);
  factory.change_flags(want_nice_coset_system ? GFF_RWS|GFF_ACCESSIBLE :
                                                 GFF_RWS|GFF_BFS|GFF_ACCESSIBLE,
                       want_nice_coset_system ? GFF_BFS : 0);
  Node_ID nid = 0;
  Node_Reference s;
  // Ensure failure state is state 0
  factory.find_state((unsigned char *) &nid,sizeof(Node_ID));

  Node_Iterator ni(nm);
  ni.begin();
  /* Now find exactly which states we will need, so that we can create
     the object as efficiently as possible. On entry we can be sure
     that NF_NEEDED flag is not set anywhere because it has the same
     value as NF_NEW, which is always set when a node is first created but
     cleared by the time the node has been fully constructed. */
  bool inside = true;
  Node_Count done = 0;
  Node_Count to_do = rm->node_count(true);
  while (ni.scan(&s,inside))
  {
    nid = s;
    inside = s->last_letter() < nr_generators;
    if (s->is_final())
    {
      if (flags & RWSC_CRITICAL_ONLY)
      {
        if (!s->flagged(EQ_INTERESTING+EQ_HAS_DIFFERENCES+
                           EQ_HAS_PRIMARY_DIFFERENCES+EQ_AXIOM))
          nid = 0;
      }
      else if (flags & RWSC_MINIMAL && !s->fast_is_primary())
        nid = 0;
      if (nid)
      {
        s = s->parent(nm);
        while (!s.is_null() && !s->flagged(NF_NEEDED))
        {
          s->node(nm).set_flags(NF_NEEDED);
          s = s->parent(nm);
        }
      }
    }
    if (!(char) done++)
      container.status(2,1,"Examining nodes (" FMT_NC " of " FMT_NC ")\n",done,to_do);
  }
  to_do = done;

  // Ensure that first equation is number 1
  nid = 0;
  eqn.find_entry((unsigned char *) &nid,sizeof(Node_ID));
  // Create the correct initial state
  ni.begin();
  ni.scan(&s);// This gets us the root */
  s->node(nm).set_flags(NF_NEEDED); // needed for when we have an RWS with no equations at all.
  nid = s;
  factory.find_state((unsigned char *) &nid,sizeof(Node_ID));

  State_ID state = 0;
  State_Count count = 2;
  State_ID ceiling = 1;
  Word_Length state_length = 0;
  Element_Count eqn_count = 1;
  const char * message;
  if (flags & RWSC_NEED_FSA)
    if (flags & RWSC_FSA_ONLY)
      message = "Building word-acceptor (" FMT_ID " of " FMT_ID " nodes to do)\n";
    else
      message = "Building index automaton (" FMT_ID " of " FMT_ID
                " nodes to do). " FMT_ID " equations\n";
  else
    message = "Sorting equations (" FMT_ID " of " FMT_ID
              " nodes to do). " FMT_ID " equations\n";
  while (factory.get_state_key(&nid,++state))
  {
    s = Node_Reference(nm,nid);
    if (!(char) state)
      container.status(2,1,message,count - state,count,eqn_count);
    for (i = 0; i < nr_generators;i++)
    {
      Ordinal g = transition_order[i];
      Node_Reference ns = nm.new_state(s,g);
      nid = ns;
      if (ns->is_final())
      {
        if (flags & RWSC_CRITICAL_ONLY)
        {
          if (!ns->flagged(EQ_INTERESTING+EQ_HAS_DIFFERENCES+
                           EQ_HAS_PRIMARY_DIFFERENCES+EQ_AXIOM))
            nid = 0;
        }
        else if (flags & RWSC_MINIMAL)
        {
          ns = ns->primary(nm,ns);
          nid = ns;
        }
        Element_ID eqn_nr = flags & RWSC_FSA_ONLY ? 0 : eqn.find_entry((unsigned char *) &nid,sizeof(Node_ID));
        if (eqn_nr >= eqn_count)
        {
          Simple_Equation se(nm,ns);
          const Word & rhs = se.get_rhs();
          const Word & lhs = se.get_lhs();
          Extra_Packed_Word packed_lhs(lhs);
          Extra_Packed_Word packed_rhs(rhs);
          packed_lhses[eqn_nr] = packed_lhs;
          packed_rhses[eqn_nr] = packed_rhs;
          eqn_count++;
        }
        transition[g] = - (State_ID) eqn_nr;
      }
      else
      {
        if (want_nice_coset_system && state != 1)
        {
          /* First test ensures that after G generators or coset symbol
             only G generators occur.
             Second test ensures that after H generators only H generators
             or coset symbol occur. So only words are G words, H words,
             or h*_H*g words */
          if (s->last_letter() <= pd.coset_symbol && g >= pd.coset_symbol ||
              s->last_letter() > pd.coset_symbol && g < pd.coset_symbol)
          {
            transition[g] = 0;
            continue;
          }
        }
        if (!ns->flagged(NF_NEEDED))
        {
          Ordinal_Word ow(rm->alphabet());
          ns->read_word(&ow,nm);
          Word_Length start = 0;
          Word_Length l = ow.length();
          do
            ns = nm.get_state(ow.buffer(),l,++start);
          while (!ns->flagged(NF_NEEDED));
          nid = ns;
        }
        transition[g] = factory.find_state((unsigned char *) &nid,sizeof(Node_ID));
        if (transition[g] >= count)
        {
          if (state >= ceiling)
          {
            state_length++;
            ceiling = transition[g];
          }
          count++;
        }
      }
    }
    if (flags & RWSC_NEED_FSA)
      factory.set_transitions(state,transition);
  }

  // now clear NF_NEEDED again, since it has the same value as NF_NEW
  ni.begin();
  inside = true;
  done = 0;
  while (ni.scan(&s,inside))
  {
    if (!s->is_final())
      s->node(nm).clear_flags(NF_NEEDED);
    inside = s->last_letter() < nr_generators;
    if (!(char) done++)
      container.status(2,1,"Examining nodes (" FMT_NC " of " FMT_NC ")\n",done,to_do);
  }

  delete [] transition_order;
  delete [] transition;
  eqn.clean(); // we need to keep this - it contains the equations
  factory.remove_keys();
  factory.change_flags(GFF_TRIM); /* we don't set the minimised flag,
                                     because although we are minimised as an
                                     RWS we are not minimised as an FSA since
                                     the accepted language is described by
                                     a word-acceptor which usually has fewer
                                     states */
  if (flags & RWSC_NEED_FSA)
  {
    if (taker)
      *taker = FSA_Factory::copy(factory);
    else
      fsa__ = FSA_Factory::copy(factory);
  }
  else
    fsa__ = 0;


  if (flags & RWSC_NEED_FSA)
    container.progress(1,"Index automaton has " FMT_ID " states and " FMT_ID
                         " equations\n",
                       factory.state_count(),eqn.count()-1);
}

/**/

Rewriting_System::Rewriting_System(const Rewriting_System & other) :
  Delegated_FSA(other.container,other.base_alphabet,true),
  confluent(other.confluent),
  inverses(0),
  source(other.source),
  eqn(* new Hash(other.equation_count(),sizeof(State)))
{
  /* This constructor should only be used to copy an RWS with an FSA
     otherwise the inverses field cannot be set */
  if (other.fsa())
  {
    fsa__ = FSA_Factory::copy(* other.fsa());
  }
  Ordinal nr_generators = base_alphabet.letter_count();
  inverses = new Ordinal[nr_generators];
  memcpy(inverses,other.inverses,nr_generators*sizeof(Ordinal));
  Element_Count count = other.equation_count();
  State s = 0;
  Ordinal_Word lhs_word(base_alphabet);
  Ordinal_Word rhs_word(base_alphabet);
  eqn.manage(packed_lhses);
  eqn.manage(packed_rhses);
  eqn.find_entry((unsigned char *) &s,sizeof(State));
  for (Element_ID eqn_nr = 1; eqn_nr < count;eqn_nr++,s++)
  {
    other.read_lhs(&lhs_word,eqn_nr);
    other.read_rhs(&rhs_word,eqn_nr);
    eqn.find_entry((unsigned char *) &s,sizeof(State));
    Extra_Packed_Word packed_lhs(lhs_word);
    Extra_Packed_Word packed_rhs(rhs_word);
    packed_lhses[eqn_nr] = packed_lhs;
    packed_rhses[eqn_nr] = packed_rhs;
  }
  eqn.clean();
}

/**/

Rewriting_System::Rewriting_System(const MAF & maf,String fsa_filename) :
  Delegated_FSA(maf.container,maf.alphabet,true),
  eqn(* new Hash(maf.axiom_count()+1,sizeof(Linked_Packed_Equation *))),
  source(0)
{
  /* this function builds a KBMAG style rewriting system from a MAF
     presentation */
  Ordinal nr_generators = maf.generator_count();

  /* Set up miscellaneous information needed for correct output */
  inverses = new Ordinal[nr_generators];
  for (Ordinal i = 0; i < nr_generators;i++)
    inverses[i] = maf.inverse(i);
  confluent = maf.properties().is_confluent;

  /* Now create the equation database */
  const Linked_Packed_Equation *axiom = 0;
  Element_Count eqn_count = 1;

  // Ensure that first equation is number 1
  eqn.manage(packed_lhses);
  eqn.manage(packed_rhses);
  eqn.find_entry((unsigned char *) &axiom,sizeof(axiom));
  for (axiom = maf.first_axiom();axiom;axiom = axiom->get_next())
  {
    Simple_Equation se(base_alphabet,*axiom);
    Element_ID eqn_nr = eqn.find_entry((unsigned char *) &axiom,sizeof(axiom));
    if (eqn_nr >= eqn_count)
    {
      const Word & rhs = se.get_rhs();
      const Word & lhs = se.get_lhs();
      Extra_Packed_Word packed_lhs(lhs);
      Extra_Packed_Word packed_rhs(rhs);
      packed_lhses[eqn_nr] = packed_lhs;
      packed_rhses[eqn_nr] = packed_rhs;
      eqn_count++;
    }
  }

  eqn.clean();
  if (fsa_filename)
    fsa__ = FSA_Factory::create(fsa_filename,&container);
  else
    fsa__ = 0;
}

/**/

void Rewriting_System::set_name(String name_)
{
  if (name)
    delete [] name;
  name = name_.clone();
  String_Buffer sb(String(name).length()+14);
  if (fsa__)
    fsa__->set_name(String::make_filename(&sb,"",name,".reductionFSA",false));
}

/**/

void Rewriting_System::print(Output_Stream *kb_stream,Output_Stream * rws_stream) const
{
  /* This function outputs a RWS in the same format as KBMAG
     (but with different whitespace).*/

  Ordinal_Word lhs_word(base_alphabet,1);
  Ordinal_Word rhs_word(base_alphabet,1);
  if (proof_trace && source)
    proof_trace->final_system(*this,*source,source->status().some_discarded ||
                                            source->status().some_ignored);
  container.output(kb_stream,"#Generated by MAF\n");
  container.output(kb_stream,"%s := rec\n"
                             "(\n"
                             "  isRWS := true,\n"
                             "  isConfluent := %s,\n",
                             name ? name : "_RWS",
                             confluent ? "true" : "false");
  container.output(kb_stream,"  generatorOrder := [");
  base_alphabet.print(container,kb_stream,APF_Bare);
  container.output(kb_stream,"],\n");
  base_alphabet.print_ordering(container,kb_stream);
  Ordinal nr_letters = base_alphabet.letter_count();
  if (nr_letters)
  {
    container.output(kb_stream,"  inverses :=       [");
    container.output(kb_stream,"%s",base_alphabet.glyph(inverses[0]).string());
    for (Ordinal i = 1; i < nr_letters;i++)
      container.output(kb_stream,",%s",base_alphabet.glyph(inverses[i]).string());
    container.output(kb_stream,"],\n");
  }
  container.output(kb_stream,"  equations := \n"
                             "  [\n");
  String_Buffer lsb,rsb;
  Element_Count count = eqn.count();
  for (Element_ID id = 1; id < count;id++)
  {
    read_lhs(&lhs_word,id);
    read_rhs(&rhs_word,id);
    base_alphabet.format(&lsb,lhs_word);
    base_alphabet.format(&rsb,rhs_word);
    if (id+1 < count)
      container.output(kb_stream,"    [%s,%s],\n",lsb.get().string(),rsb.get().string());
    else
      container.output(kb_stream,"    [%s,%s]\n",lsb.get().string(),rsb.get().string());
  }
  container.output(kb_stream,"  ]\n);\n");
  if (fsa__ && rws_stream)
    fsa__->print(rws_stream);
}

/**/

unsigned RWS_Reducer::reduce(Word * word,const Word & start_word,
                             unsigned flags,const FSA *)
{
  /* Convert word to reduced form by applying equations.
     Return value is number of times word gets reduced.
     In this code we cope only with reductions that do not increase
     the length of the word. If we encounter such a reduction we
     call hard_reduce() to complete the reduction. */

  unsigned retcode = 0;
  State_ID si, nsi;
  Ordinal rvalue = PADDING_SYMBOL; // Initialised to shut up the compiler
  Word_Length length = start_word.length();
  Word_Length valid_length = 0;
  Ordinal_Word rhs_word(rws.base_alphabet);
  State_ID eqn_nr = 0;
  if (!length)
    return 0;
  if (!state || max_stack < length)
  {
    if (state)
      delete [] state;
    state = new State_ID[(max_stack = length)+1];
  }

  if (word != &start_word)
  {
    word->set_length(length);
    word_copy(*word,start_word,length);
  }
  Ordinal * values = word->buffer();
  si = state[0] = rws.initial_state();
  if (flags & WR_PREFIX_ONLY)
  {
    rvalue = values[--length];
    word->set_length(length); /* we have to set this again in case there are no reductions */
  }
  while (valid_length < length)
  {
    nsi = state[valid_length+1] = rws.new_state(si,values[valid_length]);
    valid_length++;
    if (nsi < 0)
    {
      if (flags & WR_CHECK_ONLY)
        return 1;
      if (eqn_nr != -nsi)
      {
        eqn_nr = -nsi;
        rws.read_rhs(&rhs_word,eqn_nr);
      }
      Word_Length lhs_length = rws.lhs_length(eqn_nr);
      Word_Length rhs_length = rhs_word.length();
      if (rhs_length > lhs_length)
      {
        retcode += hard_reduce(word,flags,valid_length,rhs_word);
        break;
      }
      else
      {
        retcode++;
        valid_length -= lhs_length;
        memcpy(values+valid_length,rhs_word.buffer(),rhs_length*sizeof(Ordinal));
        if (rhs_length < lhs_length)
        {
          memcpy(values+valid_length+rhs_length,values+valid_length+lhs_length,
                 (length-valid_length-lhs_length)*sizeof(Ordinal));
          length -= lhs_length - rhs_length;
          word->set_length(length);
        }
        if (flags & WR_ONCE)
          break;
        si = state[valid_length];
      }
    }
    else
      si = nsi;
  }

  if (flags & WR_PREFIX_ONLY)
    word->append(rvalue);
  return retcode;
}

/**/

unsigned RWS_Reducer::hard_reduce(Word * word,unsigned flags,
                                  Word_Length valid_length,
                                  Ordinal_Word &rhs_word)
{
  /* Convert word to reduced form by applying equations.
     In this version we cater for word orders that can increase the length.
     We avoid moving the whole word up and down, as this is very expensive for
     long words.
     On entry we have just arrived at the first reduction that increases the
     word length. We move the unread part of the word to a stack.

     Initially we are at a reduction, we remove the part to be rewritten from
     the word, and put the replacement on the stack.
     Then we move letters one at a time from the stack to the word.
  */
  State_ID nsi = state[valid_length];
  State_ID eqn_nr = -nsi;
  stack_word.allocate(MAX_WORD,false);
  Ordinal * stack = stack_word.buffer();
  unsigned retcode = 0;
  Ordinal * values = word->buffer();
  Word_Length max_length = MAX_WORD;
  if (flags & WR_PREFIX_ONLY)
    max_length--;
  Word_Length stack_top = 0;
  Word_Length allocated = word->length(); // in fact more may already have been
                                          // allocated, but we cannot tell
  Total_Length total_length = allocated;
  bool reducing = true;

  for (Word_Length i = allocated; i > valid_length;)
    stack[stack_top++] = values[--i];

  for (;;)
  {
    if (nsi < 0)
    {
      if (eqn_nr != -nsi)
      {
        eqn_nr = -nsi;
        rws.read_rhs(&rhs_word,eqn_nr);
      }
      Word_Length rhs_length = rhs_word.length();
      Word_Length lhs_length = rws.lhs_length(eqn_nr);
      total_length -= lhs_length - rhs_length;
      valid_length -= lhs_length;
      if (stack_top + rhs_length > MAX_WORD)
        MAF_INTERNAL_ERROR(rws.container,
                           ("Maximum word length exceeded in"
                            " RWS_Reducer::hard_reduce()\n"));
      const Ordinal * rvalues = rhs_word.buffer();
      for (Word_Length i = rhs_length; i > 0;)
        stack[stack_top++] = rvalues[--i];
      retcode++;
      if (flags & WR_ONCE)
        reducing = false;
    }

    if (stack_top)
    {
      if (valid_length == max_stack)
      {
        State_ID * new_state = new State_ID[(max_stack = MAX_WORD)+1];
        memcpy(new_state,state,(valid_length+1)*sizeof(State_ID));
        delete [] state;
        state = new_state;
      }
      if (valid_length == max_length)
        MAF_INTERNAL_ERROR(rws.container,
                           ("Maximum word length exceeded in"
                            " RWS_Reducer::hard_reduce()\n"));
      if (valid_length >= allocated)
      {
        if (total_length > MAX_WORD)
        {
          /* we deliberately don't immediately detect situation where maximum word
             length has been exceeded, because we may recover before
             either the stack or the answer gets too long.*/
          word->allocate(allocated = MAX_WORD,true);
        }
        else
          word->allocate(allocated = Word_Length(total_length),true);
        values = word->buffer();
      }
      Ordinal g = values[valid_length] = stack[--stack_top];
      if (reducing)
        nsi = state[valid_length+1] = rws.new_state(state[valid_length],g);
      else
        nsi = 0;
      valid_length++;
    }
    else
      break;
  }

  word->set_length(valid_length);
  return retcode;
}

/**/

Rewriting_System * Rewriting_System::create(String rws_filename,
                                            String fsa_filename,
                                            Container * container,
                                            bool must_succeed)
{
  Rewriting_System * rws = 0;
  for (;;)
  {
    unsigned flags = must_succeed ? 0 : CFI_ALLOW_CREATE;
    if (fsa_filename==0)
      flags |= CFI_CREATE_RM;
    else
      flags |= CFI_RAW;
    MAF & maf = * MAF::create_from_rws(rws_filename,container,flags);
    if (!&maf)
      return 0;
    if (fsa_filename)
      rws = new Rewriting_System(maf,fsa_filename);
    else
      rws = new Rewriting_System(&maf.rewriter_machine(),false);
    delete &maf;
    if (!rws->fsa() && fsa_filename)
    {
      /* In this case we could not read the .reduce file. We can recreate the
         FSA by building a Rewriter_Machine from the axioms and then create
         a Rewriting_System from that. */
      fsa_filename = 0;
      delete rws;
      rws = 0;
    }
    else
      break;
  }
  return rws;
}
