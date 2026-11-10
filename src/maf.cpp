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


// $Log: maf.cpp $
// Revision 1.19  2010/12/09 11:16:31Z  Alun
// Now works when an unnamed RWS is saved
// Revision 1.18  2010/07/08 11:10:24Z  Alun
// save() changed to allow output to stdout with null filename
// Revision 1.17  2010/06/10 13:57:27Z  Alun
// All tabs removed again
// Revision 1.16  2010/05/17 00:03:03Z  Alun
// order() method added. Various new automata added to the GAT enumeration.
// save_fsa() now handles subgroup word acceptor better, to improve KBMAG
// compatibility
// Revision 1.15  2009/10/13 22:43:28Z  Alun
// Compatibility of FSA utilities with KBMAG improved.
// Compatibility programs produce fewer FSAs KBMAG would not (if any)
// Revision 1.14  2009/09/13 11:36:37Z  Alun
// printf() style functions cleaned up considerably
// code now compiles cleanly in g++ both 32-bit and 64-bit
// all tabs removed from source files
// Revision 1.13  2009/09/01 14:24:40Z  Alun_Williams
// Long comment updated
// Revision 1.12  2008/11/03 10:38:14Z  Alun
// Re-check in of source from web site after catastrophic loss of correct state
// of RCS archive.
// Revision 1.13  2008/11/03 11:38:13Z  Alun
// Early November checkin
// Completely reworked pool and Difference_Tracker
// Ported to Darwin
// Revision 1.12  2008/10/09 15:52:11Z  Alun
// Ensure validating word-difference machine is not missing initial g,g transitions
// Revision 1.11  2008/09/22 09:01:13Z  Alun
// Final version built using Indexer.
// Revision 1.6  2007/12/20 23:27:31Z  Alun
//

/* The user must first create a MAF object and set up the presentation of
   the object to be worked on. This will usually be done by running a
   program such as automata that reads a file containing a
   "rewriting system", but it may be done by directly calling MAF APIs.

   The object may have a previously computed structure, in which case
   various FSA describing this will have been computed, and we can simply
   load these FSA and compute with them.

   Or, more commonly, the presentation describes an unknown group or monoid
   whose structure we must deduce from the axioms. In this case the
   presentation must usually be turned into a Rewriter_Machine instance.
   Rewriter_Machine is basically a class for processing rewriting systems,
   but one in which the states and equations are labelled with pointers
   rather than integers. MAF will then use the methods provided by
   Rewriter_Machine to analyse the presentation and do one, or possibly two,
   of the following three things:

   1) It will successfully create a confluent rewriting system (explained
   below) for the group or monoid in KBMAG format. In this case we may also
   do 2) as well, unless the options supplied tell us not to.

   2) If the object described by the presentation is a group and the set of
   generators is closed under inversion, and the group is automatic with
   respect to the the chosen generators and word ordering method then the
   automatic structure is computed. This includes:
   i) various word-difference machines
   ii) the word-acceptor
   iii) the general multiplier.

   The second of these is usually the most useful, as it can be used to
   decide in linear time whether any particular word is reduced or not.
   However to reduce a word that is not accepted, the other one of the other
   machines is needed, and word reduction is quadratic.

   We may also do 2) if we are analysing a coset system.

   3) The program will run for ever until you stop it, or it has crashed
   because it has used up all the memory, or its internal limits have
   been exceeed.

   Other options
   =============

   Code to perform Tietze tranformations and do coset enumeration is
   now available, and this does not require a Rewriter_Machine (though
   the former option may well create one in a subsidiary object).

   Confluent Rewriting Systems
   ===========================

   A rewriting machine consists of a set of rewriting rules, that is
   pairs of words in which the first word, or left hand side (LHS), is to
   be replaced by the the second word (RHS), together with an automaton
   which accepts characters from a word one at a time, and can recognise
   words that are an LHS and "rewrite" the word by performing the required
   replacment. When no more rules can be applied to the input the resulting
   word is said to be reduced.

   A rewriting system generally does not give a unique possible reduction for
   every word, since if a word contains subwords for two different LHS
   then depending on which we perform first we may get a different answer.
   Such pairs of words must be equal in the group or monoid described by
   a rewriting system.  If in fact every possible path of reductions on every
   possible input word leads to the same final word the rewritig system is
   said to be confluent.

   MAF tries to create a confluent rewriting system machine corresponding to
   the object described by the axioms. Rewriter_Machine modifies the initial
   rewriting system that is created by by inserting equations corresponding
   to axioms into the the (trivial) rewriting system for the free monoid on
   the same symbols, by looking for new rules which can be deduced by finding
   two different reductions for the same input, which amounts to applying the
   associative law to this input. In order to do this the program repeatedly
   considers all "plausible pairs" of equations, until no more pairs can be
   found. If this desirable state of affairs arises the presentation is
   called "confluent".

   A pair is plausible if a suffix of the lhs of the first equation is a
   prefix of the second's lhs, and the pair has not already been considered in
   this order. (In fact the same pair of equations (in the same order) can
   be plausible in more than one way, so these "pairs" are really triples,
   consisting of the two equations and the length of the "overlap").

   As each new equation is inserted some further modifications are made to
   ensure a minimum level of consistency in the tree.

   For more details see maf_rm.cpp and the modules it refers to.

   This process is called the Knuth Bendix procedure, and it can
   take an arbitrarily long time, and may not ever terminate for infinite
   groups and monoids, because each new equation that is added creates yet
   more pairs of equations to be considered.

   Now, if a group has an shortlex automatic structure, then it can be shown
   there is a finite set of words, called the "word-differences" formed as
   follows.
   Suppose u and v are accepted words and Uv=x (with x a generator or the
   identity). The set of word-differences is the set of all group elements
   of the form p(u)^-1*p(v) where p(u) and p(v) are prefixes of the same
   length of u and v (or the whole of the shorter word in the case where
   one of the words is shorter than the other and the prefix length of the
   longer word is already longer than the whole of the shorter word). If
   we know this set we can use it to calculate the automata described above
   for the group, and when we have done so can prove that the automata
   satisfy the group axioms.
   If we know only a subset of this set, or are accepting two different
   words for the same group element in our difference machine, then when we
   try to perform this calculation something will go wrong. Also, we will
   usually end up using much more memory than we would have if we had started
   from the correct set of differences. The word-acceptor from such
   a failed calculation is of some use however - it always accepts at least
   one word for every element in the group

   The program works roughly as follows:

   0  Before each pass of KB completion MAF will decide whether to try
   to build the automata from its current best guess at the word
   differences. Usually it will only do so when:
   i) the size of this set appears to have more or less stabilised and
   ii) the new equations being formed in KB expansion are long in
   comparison to the longest known word-difference and
   iii) the "height" of the Rewriter_Machine is still increasing.

   If i) or ii) have not happened yet then we probably do not yet know
   all the word-differences. If iii) is not true, then it is quite likely that
   the presentation is confluent, so it is worth continuing the KB process.

   This step is only performed when the presentation specified a set of
   generators closed under inversion, and the option to iterate until
   confluence was not set.

   1  All the equations are kept in a tree that effectively acts as an FSA,
   and every time an equation is added, any simplifications of existing
   equations that are possible are discovered immediately and performed as
   already outlined.

   2 The process of discovering pairs of equations for KB also uses the tree
   structure to avoid the need for explicit string comparison. This saves
   a great deal of time when there are very many equations and not too many
   generators. In any one pass this program does not consider all possible
   pairs, but only pairs which satisfy certain criteria. These criteria are
   gradually relaxed, so that eventually, in the confluent case, all plausible
   pairs have been considered. In the case where confluence is never going to
   be achieved we hope our criteria have helped us find the word-differences
   in the minimum amount of time.

   3 The program also forms equations in which no prefix of the LHS is
   reducible and the lhs and rhs have no common prefix or suffix.  Such
   equations are not needed for the minimal rewriting system, but are needed
   for calculating the word-differences. Such equations do not need to be
   considered as the left partner of a plausible pair. However, they may
   be considered on the right in this implementation, though this depends
   on configuration options, and would only be useful for non-geodesic
   word orderings.

   4 When new equations are first created, unless they satisfy a length
   criterion, they are not considered when creating new equations (apart from
   in the simplification process), and are only "waiting" equations. Each
   pass of the algorithm first creates a candidate list of equations to be
   considered, then considers this list in its entirety (possibly adding a few
   new equations along the way if they are short enough), and then finally
   inspects the current list of pending equations, and adds certain of them
   (usually the shortest) to the list of "approved" equations that will
   be considered in the next KB pass. There are two classes of such equation:
   equations which do not increase the height of the tree by more than one
   are inserted into the tree at once, equations longer than this are generally
   placed in a "pool" and are only put into the tree when it has grown tall
   enough for them to be worth considering (though some longer equations will
   often be inserted, depending on what options are in effect)

   No equation is ever completely discarded until it has been proved to have a
   reducible prefix, and the reduction of the LHS and RHS do not add to our
   knowledge. This means the results of the program are definitely correct. Or,
   more accurately, MAF knows when the equation arising from an overlap might
   have been discarded, and will recreate it on a later pass.
*/

#include "mafauto.h"
#include "mafword.h"
#include "container.h"
#include "maf_rm.h"
#include "maf_rws.h"
#include "maf_dr.h"
#include "maf_we.h"
#include "rubik.h"
#include "mafcoset.h"
#include "maf_trace.h"

bool Word_Reducer::reducible(const Word & word,Word_Length length)
{
  Ordinal_Word test(Subword((Word &) word,0,length));
  return reduce(&test,test,WR_CHECK_ONLY)!=0;
}

Container * MAF::create_container(Platform * platform)
{
  return Container::create(platform);
}

MAF * MAF::create(Container * container,
                  const Options * options,
                  Alphabet_Type alphabet_type,
                  Presentation_Type presentation_type)
{
  bool auto_delete(container==0);
  if (!container)
    container = create_container();
  return new MAF(*container,alphabet_type,presentation_type,auto_delete,options);
}

/**/

MAF::MAF(Container & container,Alphabet_Type alphabet_type,
         Presentation_Type presentation_type,bool delete_container,
         const Options * options_) :
  Presentation(container,alphabet_type,presentation_type,delete_container),
  aborting(false),
  fsas(real_fsas),
  rm(0),
  automata(0),
  wr(0),
  provisional_wr(0),
  validator(0)
{
  if (options_)
    options = *options_;
}

/**/

MAF::~MAF()
{
  if (automata)
    delete automata;
  if (rm)
    delete rm;
  if (wr)
    delete wr;
  if (provisional_wr)
    delete provisional_wr;
  if (validator)
  {
    delete (FSA *) &validator->dm;
    delete validator;
    validator = 0;
  }
}

/**/

// methods for processing words

bool MAF::reduce(String_Buffer * rword,String word) const
{
  bool retcode = false;
  Ordinal_Word *parsed = parse(word);
  if (parsed)
  {
    reduce(parsed,*parsed);
    parsed->format(rword);
    delete parsed;
    retcode = true;
  }
  return retcode;
}

/**/

bool MAF::reducible(String word,String_Length length) const
{
  bool retcode = false;
  Ordinal_Word *parsed = parse(word,length);
  if (parsed)
  {
    retcode = reducible(*parsed,parsed->length());
    delete parsed;
  }
  return retcode;
}

/**/

bool MAF::L2_accepted(String word,String_Length length) const
{
  bool retcode = false;
  Ordinal_Word *parsed = parse(word,length);
  if (parsed)
  {
    retcode = L2_accepted(*parsed);
    delete parsed;
  }
  return retcode;
}

/**/

bool MAF::L2_accepted(const Word &word) const
{
  return !reducible(word,word.length()-1);
}

/**/

bool MAF::insert_axiom(const Word & lhs,const Word & rhs,unsigned flags)
{
  if (!rm)
    realise_rm();

  if (rm->add_axiom(lhs,rhs,flags))
  {
    Ordinal_Word lhs_word(lhs);
    Ordinal_Word rhs_word(rhs);

    for (Generator_Permutation * p = permutation;p;p = p->next)
    {
      Word_Length i,j;
      const Ordinal_Word & perm = *p->word;
      Word_Length pl = perm.length();
      Word_Length l = lhs_word.length();

      for (i = 0;i < l;i++)
        for (j = 0;j < pl;j += 2)
          if (perm.value(j) == lhs_word.value(i))
          {
            lhs_word.set_code(i,perm.value(j+1));
            break;
          }
      l = rhs_word.length();
      for (i = 0;i < l;i++)
        for (j = 0;j < pl;j += 2)
          if (perm.value(j) == rhs_word.value(i))
          {
            rhs_word.set_code(i,perm.value(j+1));
            break;
          }
      insert_axiom(lhs_word,rhs_word,0);
    }
    return true;
  }

  return false;
}

/**/

unsigned MAF::reduce(Word * rword,const Word & word,unsigned flags) const
{
  if (wr)
    return wr->reduce(rword,word,flags);
  if (rm)
  {
    if (rm->status().want_differences && !provisional_wr)
      provisional_wr = new Strong_Diff_Reduce(rm);
    if (provisional_wr)
      return provisional_wr->reduce(rword,word,flags);
    return rm->reduce(rword,word,flags);
  }

  return 0;
}

/**/

bool MAF::reducible(const Word &word,Word_Length length) const
{
  if (wr)
    return wr->reducible(word,length);
  if (automata)
  {
    if (length == WHOLE_WORD)
      length = word.length();
    return automata->word_accepted(word.buffer(),length);
  }
  if (rm)
  {
    if (rm->status().want_differences && !provisional_wr)
      provisional_wr = new Strong_Diff_Reduce(rm);
  }
  if (provisional_wr)
    return provisional_wr->reducible(word,length);
  return 0;
}

/**/

void MAF::realise_rm()
{
  if (!rm)
  {
    if (!flags_set)
      set_flags();
    if (proof_trace && !proof_trace->started())
    {
      /* The header has to come before the first axiom */
      if (!proof_trace->start((MAF &) *this))
        exit(1);
    }
    rm = new Rewriter_Machine((MAF &)*this);
    /* In this case any axioms added so far have been stored in the
       presentation */
    for (Ordinal g = 0; g < nr_generators;g++)
      if (inverse(g) != INVALID_SYMBOL)
      {
        Ordinal_Word s(alphabet,2);
        Ordinal_Word id_word(alphabet);
        s.set_code(0,g);
        s.set_code(1,inverse(g));
        insert_axiom(s,id_word,AA_ADD_TO_RM);
      }
    for (const Linked_Packed_Equation * axiom = Presentation::first_axiom(false);axiom;axiom = axiom->get_next())
    {
      Simple_Equation se(alphabet,*axiom);
      add_axiom(se.lhs_word,se.rhs_word);
    }
  }
}

/**/

void MAF::grow_automata(FSA_Buffer * buffer,
                        unsigned save_flags,
                        unsigned retain_flags,
                        unsigned exclude_flags)
{
  if (!automata)
  {
    realise_rm();
    rm->start();
    int action;
    automata = new Group_Automata;
    while ((action = rm->expand_machine())!=0)
    {
      if (automata->build_vital(rm,false,action))
      {
        if (aborting || action != 2)
          break;
        /* In this case we are periodically saving where we got to so far */
        automata->grow_automata(rm,buffer,GA_PDIFF2|GA_DIFF2|GA_RWS,0,exclude_flags);
      }
      automata->erase();
    }
    if (action == 0)
      automata->build_vital(rm,true,action);
  }
  automata->grow_automata(rm,buffer,save_flags,retain_flags,exclude_flags);
  automata->transfer(&real_fsas);
}

/**/

void MAF::set_validation_fsa(FSA_Simple * fsa)
{
  if (validator)
  {
    delete (FSA *) &validator->dm;
    delete validator;
  }
  Ordinal nr_symbols = fsa->base_alphabet.letter_count();
  for (Ordinal g = 0;g < nr_symbols;g++)
    fsa->set_transition(1,fsa->base_alphabet.product_id(g,g),1);
  validator = new Diff_Reduce(fsa);
}

/**/

bool MAF::is_valid_equation(String *reason,const Word &lhs,const Word & rhs) const
{
  if (is_rubik)
  {
    String_Buffer sb1,sb2;
    lhs.format(&sb1);
    rhs.format(&sb2);
    if (!Rubik::compare(sb1.get(),sb2.get()))
    {
      *reason = "Untrue Rubik equation\n";
      return false;
    }
  }

  if (validator)
  {
    Ordinal_Word lhs_word(group_alphabet());
    Ordinal_Word rhs_word(group_alphabet());
    group_word(&lhs_word,lhs);
    group_word(&rhs_word,rhs);
    if (!validator->dm.product_accepted(lhs_word,rhs_word))
    {
      validator->reduce(&lhs_word,lhs_word);
      validator->reduce(&rhs_word,rhs_word);
      if (!validator->dm.product_accepted(lhs_word,rhs_word))
      {
        *reason = "Equation is not accepted by validating word-difference machine\n";
        return false;
      }
    }
  }
  return true;
}

/**/

void MAF::save(String save_filename)
{
  if (filename != save_filename && save_filename)
    filename = save_filename.clone();
  /* If no name has been set yet we had better create one,
     because otherwise any FSA we save will have an invalid name */
  if (!name)
    name = String(is_coset_system ? "_RWS_Cos" : "_RWS").clone();
  Output_Stream * os = container.open_text_output_file(save_filename ? filename : save_filename);
  print(os);
  container.close_output_file(os);
}

/**/

void MAF::print(Output_Stream * stream) const
{
  if (rm)
  {
    Rewriting_System rws(rm,RWSC_MINIMAL);
    if (name != 0)
      rws.set_name(name);
    rws.print(stream,(Output_Stream *) 0);
  }
  else
  {
    Rewriting_System rws(*this);
    if (name != 0)
      rws.set_name(name);
    rws.print(stream,(Output_Stream *) 0);
  }
}

/**/

bool MAF::output_gap_presentation(Output_Stream * stream,bool change_alphabet) const
{
  /* This converts our presentation into something GAP can understand.
     This method is primarily intended for KBMAG compatibility, so that
     we have an easy way to output subgroup presentations in the GAP
     format it uses.
     The change_alphabet flag makes the generators look like _x1,_x2,....
     This is not very nice, so we try to preserve as many generator names
     as we can if this flag is false.
     We assume that out of each pair of inverse generators at least one has a
     name not ending in ^-1. We pick this one, or if neither end in ^-1 we
     pick the first, which is called x say. The other generator is then named
     x^-1.
  */

  if (!is_group)
    return false;
  Alphabet & new_alphabet = *Alphabet::create(AT_String,container);
  Ordinal_Word relator(new_alphabet,2);
  bool started = false;
  unsigned nr_free_generators = 0;

  if (change_alphabet)
  {
    new_alphabet.set_nr_letters(nr_generators);
    new_alphabet.set_similar_letters("_x",nr_free_generators = nr_generators);
  }
  else
  {
    String_Buffer sb;

    for (Ordinal g = 0;g < nr_generators;g++)
    {
      if (g)
        sb.append(",");
      Ordinal ig = inverse(g);
      if (g <= ig)
      {
        Glyph s = alphabet.glyph(g);
        size_t l = s.length();
        if (l < 3 || strcmp(s+l-3,"^-1")!=0)
        {
          nr_free_generators++;
          sb.append(s);
        }
        else
        {
          sb.append(alphabet.glyph(ig));
          sb.append("^-1");
        }
      }
      else
      {
        Glyph s = alphabet.glyph(ig);
        size_t l = s.length();
        if (l > 3 && strcmp(s+l-3,"^-1")==0)
        {
          nr_free_generators++;
          sb.append(alphabet.glyph(g));
        }
        else
        {
          sb.append(s);
          sb.append("^-1");
        }
      }
    }
    new_alphabet.set_letters(sb.get());
  }

  container.output(stream,"_xg := FreeGroup(%d);\n",nr_free_generators);
  int j = 1;
  for (Ordinal g = 0;g < nr_generators;g++)
  {
    Glyph s = new_alphabet.glyph(g);
    size_t l = s.length();
    if (l < 3 || strcmp(s+l-3,"^-1"))
      container.output(stream,"%s := _xg.%d;\n",new_alphabet.glyph(g).string(),j++);
  }
  container.output(stream,"_x_relators :=\n[\n");

  String_Buffer sb;
  for (Ordinal g = 0;g < nr_generators;g++)
  {
    if (started)
      container.output(stream,",\n");
    Ordinal ig = inverse(g);
    if (change_alphabet || g == ig)
    {
      relator.set_code(0,g);
      relator.set_code(1,ig);
      relator.format(&sb);
      container.output(stream,"  %s",sb.get().string());
      started = true;
    }
  }

  for (const Linked_Packed_Equation * axiom = first_axiom();axiom;
       axiom = axiom->get_next())
  {
    Simple_Equation se(new_alphabet,*axiom);
    if (se.relator(&relator,*this))
    {
      relator.format(&sb);
      if (started)
        container.output(stream,",\n");
      container.output(stream,"  %s",sb.get().string());
    }
    started = true;
  }
  container.output(stream,started ? "\n];\n" : "];\n" );
  delete &new_alphabet;
  return true;
}

/**/

void MAF::save_fsa(FSA * fsa,Group_Automaton_Type gat,unsigned format_flags) const
{
  String suffix;
  String filename = this->filename;
  String name = this->name;
  String_Buffer sb;
  if (is_coset_system)
  {
    suffix = FSA_Buffer::coset_suffixes[gat];
    if (gat == GAT_Subgroup_Word_Acceptor)
    {
      filename = this->subgroup_filename;
      name = "_RWS.wa";
    }
  }
  else
    suffix = FSA_Buffer::suffixes[gat];
  fsa->save_as(filename,name,suffix,format_flags);
}

/**/

void MAF::import_difference_machine(const FSA_Simple & dm)
{
  rm->start();
  rm->import_difference_machine(dm);
}

/**/

FSA_Simple *MAF::translate_acceptor(const MAF & maf_start,
                                    const FSA & fsa_start,
                                    const Word_List &xlat_wl) const
{
  return translate_acceptor(*this,maf_start,fsa_start,xlat_wl);
}

/**/

bool MAF::polish_equation(Simple_Equation * se) const
{
  Working_Equation we(*rm,se->get_lhs(),se->get_rhs(),Derivation(BDT_Unspecified));
  if (!we.normalise(AE_KEEP_LHS))
    return false;
  se->lhs_word = we.get_lhs();
  se->rhs_word = we.get_rhs();
  reduce(&se->lhs_word,se->lhs_word,WR_PREFIX_ONLY);
  reduce(&se->rhs_word,se->rhs_word);
  return se->lhs_word != se->rhs_word;
}

/**/

bool MAF::load_reduction_method(Group_Automaton_Type flag)
{
  const FSA * fsa = 0;
  if (wr)
  {
    delete wr;
    wr = 0;
  }

  switch (flag)
  {
    case GAT_Auto_Select:
    {
      const Group_Automaton_Type path[] =
      {
        GAT_Coset_Table,
        GAT_Fast_RWS,
        GAT_Minimal_RWS,
        GAT_General_Multiplier,
        GAT_Correct_Difference_Machine,
        GAT_Full_Difference_Machine,
        GAT_Primary_Difference_Machine,
        GAT_Provisional_DM2,
        GAT_Provisional_DM1,
        GAT_Provisional_RWS,
        GAT_Auto_Select
      };

      for (int i = 0; path[i] != GAT_Auto_Select;i++)
        if (load_reduction_method(path[i]))
          return true;
      return false;
    }

    case GAT_Minimal_RWS:
    case GAT_Fast_RWS:
    case GAT_Provisional_RWS:
      fsa = load_fsas(1 << flag);
      if (fsa)
        wr = new RWS_Reducer(* (Rewriting_System *) fsa);
      return wr != 0;

    case GAT_Coset_Table:
      fsa = load_fsas(1 << flag);
      if (fsa)
        wr = new Coset_Table_Reducer(*this,*fsa);
      return wr != 0;

    case GAT_Full_Difference_Machine:
    case GAT_Primary_Difference_Machine:
    case GAT_Correct_Difference_Machine:
    case GAT_Provisional_DM1:
    case GAT_Provisional_DM2:
      fsa = load_fsas(1 << flag);
      if (fsa)
        wr = new Diff_Reduce(fsa);
      return wr != 0;

    case GAT_General_Multiplier:
      if (load_fsas(1 << flag))
      {
        fsas.gm->ensure_dense();
        wr = new GM_Reducer(*fsas.gm);
      }
      return wr != 0;

    case GAT_Deterministic_General_Multiplier:
    case GAT_GM2:
    case GAT_DGM2:
    case GAT_WA:
    case GAT_L1_Acceptor:
    case GAT_Primary_Recogniser:  // In principle we could use the three
    case GAT_Equation_Recogniser: // recognisers with Diff_Reduce
    case GAT_Reduction_Recogniser:
    case GAT_Subgroup_Word_Acceptor:
    case GAT_Subgroup_Presentation:
    case GAT_Conjugacy_Class:
    case GAT_Conjugator:
    default:
      MAF_INTERNAL_ERROR(container,
                         ("Invalid flag %d passed to "
                          "MAF::load_reduction_method()\n",flag));
      return false;
  }
}

/**/

String MAF::multiplier_name(String_Buffer * sb,const Word & word,bool coset_multiplier) const
{
  sb->set(coset_multiplier ? "mim" : "m");
  String_Buffer temp;
  Word_Length l = word.length();
  const Ordinal * values = word.buffer();
  if (l)
  {
    for (Word_Length i = 0; i < l;i++)
    {
      temp.format(i == 0 ? "%d" : "_%d",values[i]+1);
      sb->append(temp.get());
    }
  }
  else
    sb->append("0");
  return sb->get();
}

/**/

unsigned long MAF::order(Word *test,const FSA * wa) const
{
  /* This method will eventually succeed if test is a torsion element,
     but if not it may not be able to spot this. For example we
     can't even spot that a word a*b has infinite order when <a,b> is
     free abelian.
     In this kind of situation we might need to create a presentation
     for the subgroup generated by the element and show that this in
     infinite cyclic.
  */

  unsigned long answer = 1;
  reduce(test,*test);

  if (!test->length())
    return answer;

  Ordinal_Word prefix(*test);
  Ordinal_Word new_test(alphabet);
  State_ID si = wa ? wa->initial_state() : 0;

  for (;;)
  {
    /* Cyclically reduce prefix, because otherwise we have no chance of
       showing it has infinite order. */
    Ordinal c;
    Word_Length l = prefix.length();

    if (wa)
    {
      while (l > 1 && prefix.value(0) == inverse(c = prefix.value(l-1)))
      {
        new_test.set_length(0);
        new_test.append(c);
        new_test += *test;
        new_test.append(prefix.value(0));
        reduce(test,new_test);
        prefix = Subword(prefix,1,l-1);
        l -= 2;
      }

      if (wa->is_repetend(si,prefix))
      {
        answer = 0;
        return answer;
      }
    }

    for (;;)
    {
      answer++;
      prefix += *test;
      if (reduce(&prefix,prefix))
      {
        if (!prefix.length())
          return answer;
        if (wa)
          break;
      }
    }

    /* we may have test = W x w where x < test but test < W x w
      Then for some n we may get test^n = W x^n w with W x^n w now less
      than x^n, which suggests that we should try to find the order of x
      instead. However, x may in turn turn out to be V y v. In fact
      it is perfectly possible for there to be a cycle of such relations.

      Suppose we start off with bbb and this is equal to Acca,
      and that also cca=bcc
      Then we get bbbbbb=AccaAcca = Acccca = Accbcc and this is less,
      and looks cyclically reduced.
      We replace prefix by ccbccA=cccc and word becomes abbbA = cc.
      This shows why changing prefix might improve our chances of getting
      a torsion free element.
    */
    c = prefix.value(0);
    if (inverse(c) != INVALID_SYMBOL)
    {
      new_test.set_length(0);
      new_test.append(inverse(c));
      new_test += *test;
      new_test.append(c);
      reduce(test,new_test);
      prefix = Subword(prefix,1);
      prefix.append(c);
      reduce(&prefix,prefix);
      if (container.status_needed(1))
      {
        String_Buffer sb1,sb2;
        test->format(&sb1);
        prefix.format(&sb2);
        container.progress(2,"Now testing w=%s\nw^%lu=%s\n",
                             sb1.get().string(),answer,sb2.get().string());
      }
    }
  }
}
