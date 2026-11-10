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
$Log: maf.h $
Revision 1.23  2011/06/11 12:03:46Z  Alun
added new translate_acceptor() method which does things properly
Revision 1.22  2010/06/17 10:08:12Z  Alun
Some options fields were not initialised
Revision 1.21  2010/06/10 13:58:01Z  Alun
All tabs removed again
Revision 1.20  2010/06/08 06:21:28Z  Alun
Jun 2010 version. Many new/changed methods
Revision 1.19  2009/11/10 08:31:43Z  Alun
Various new options added
Revision 1.18  2009/10/13 20:46:33Z  Alun
Compatibility of FSA utilities with KBMAG improved.
Compatibility programs produce fewer FSAs KBMAG would not (if any)
Revision 1.17  2009/09/13 09:05:05Z  Alun
printf() style functions cleaned up considerably
code now compiles cleanly in g++ both 32-bit and 64-bit
all tabs removed from source files
Revision 1.16  2009/09/01 12:08:11Z  Alun_Williams
Comment added
Revision 1.15  2008/12/27 11:21:46Z  Alun
Methods to support abelianisation added
Revision 1.15  2008/10/22 01:41:12Z  Alun
Early November checkin
Completely reworked pool and Difference_Tracker
Ported to Darwin
Revision 1.14  2008/10/08 00:10:48Z  Alun
force_multiplier option added. rewriter_machine() invokes realise_rm()
Revision 1.13  2008/09/25 09:28:01Z  Alun
Final version built using Indexer.
Revision 1.7  2007/12/20 23:25:45Z  Alun
*/
#pragma once
#ifndef MAF_INCLUDED
#define MAF_INCLUDED 1
#ifndef AWCC_INCLUDED
#include "awcc.h"
#endif

/* This header file contains the public C++ interface for MAF (Monoid Automata
Factory) */

// classes declared in this header
class MAF;
class Presentation;
class Presentation_Data;
class Multiplier;
class General_Multiplier;
class GM_Reducer;

// classes referred to but defined elsewhere
class Alphabet;
class Container;
class Diff_Reduce;
class FSA;
class FSA_Simple;
class FSA_Common;
class General_Multiplier;
struct Generator_Permutation;
struct TC_Enumeration_Options;
class Group_Automata;
class Linked_Packed_Equation;
class Node_Manager;
class Ordinal_Word;
class Packed_Equation_List;
class Packed_Word;
class Platform;
class Rewriter_Machine;
class Rewriting_System;
class Simple_Equation;
class Subalgebra_Descriptor;
class Word;
class Word_Collection;
class Word_List;
class Parse_Error_Handler;
class Sorted_Word_List;

#ifndef MAFBASE_INCLUDED
#include "mafbase.h"
#endif

/*****************************************************************************/
/* CONSTANTS */

/* Codes used to represent various special subsets of words in an alphabet */
/* L0 is the language of reduced words */
const Language language_L0 = 0;
/* L1 is the language of reduced words and words that have no reducible
   subwords. This is the language of the minimal rewriting system's left
   hand sides. */
const Language language_L1 = 1;
/* L2 is the language L1 plus the words that have no reducible prefix and in
   which the only reducible subword is the suffix.
*/
const Language language_L2 = 2;
/* L3 is the language L1 plus the words that have no reducible prefix, but
  in which any suffix may be reducible. This is usually much larger than L2.
*/
const Language language_L3 = 3;
/* A is the language of all words in the alphabet */
const Language language_A = 4;

/* The Presentation_Type value in Presentation_Data determines what aspects
   of processing are carried out. Where input is read from a KBMAG format
   file generally MAF can infer the correct type to use, and you can use
   PT_General to allow this to happen. This is also true for presentations
   created directly with APIs. However, if you want to analyze a monoid with
   cancellation (e.g. the +ve integers under multiplication) you must specify
   this explicitly. Similarly you must specify PT_Rubik if you want to use the
   built in Rubik cube simulation.

   Where the details of the presentation do not conform to the flag specified
   MAF will issue warnings and change the type to match the actual
   presentation. */

enum Presentation_Type
{
  PT_General,
  PT_Monoid_Without_Cancellation,
  PT_Monoid_With_Cancellation,
  PT_Group,
  PT_Rubik,
  PT_Simple_Coset_System,
  PT_Coset_System_With_Generators,
  PT_Coset_System_With_Inverses
};

enum Compositor_Algorithm
{
  CA_Serial,
  CA_Hybrid,
  CA_Parallel
};

/* Flags to pass to grow_automata() to request various
   automata relevant to a MAF presentation */
const unsigned GA_ALL = ~0u;
const unsigned GA_RWS = 1 << GAT_Minimal_RWS;
const unsigned GA_MAXRWS = 1 << GAT_Fast_RWS;
const unsigned GA_WA = 1 << GAT_WA;
const unsigned GA_COSETS = 1 << GAT_Coset_Table;
const unsigned GA_CCLASS = 1 << GAT_Conjugacy_Class;
const unsigned GA_CONJUGATOR = 1 << GAT_Conjugator;
const unsigned GA_DIFF2 = 1 << GAT_Full_Difference_Machine;
const unsigned GA_GM = 1 << GAT_General_Multiplier;
const unsigned GA_DGM = 1 << GAT_Deterministic_General_Multiplier;
const unsigned GA_GM2 = 1 << GAT_GM2;
const unsigned GA_DGM2 = 1 << GAT_DGM2;
const unsigned GA_MINRED = 1 << GAT_L1_Acceptor;
const unsigned GA_L1_ACCEPTOR = GA_MINRED;
const unsigned GA_MINKB = 1 << GAT_Primary_Recogniser;
const unsigned GA_DIFF1C = 1 << GAT_Primary_Difference_Machine;
const unsigned GA_MAXKB = 1 << GAT_Equation_Recogniser;
const unsigned GA_RR = 1 << GAT_Reduction_Recogniser;
const unsigned GA_DIFF2C = 1 << GAT_Correct_Difference_Machine;
const unsigned GA_SUBWA = 1 << GAT_Subgroup_Word_Acceptor;
const unsigned GA_SUBPRES = 1 << GAT_Subgroup_Presentation;
const unsigned GA_PDIFF2 = 1 << GAT_Provisional_DM2;

class FSA_Buffer
{
  public:
    static String const rws_suffixes[];
    static String const suffixes[];
    static String const coset_suffixes[];
    /* The member variables must be laid out in the order of the
       GAT_xx enumeration. If this is changed then mafload.cpp
       will need to be changed not to rely on this */
    Rewriting_System * min_rws;
    Rewriting_System * fast_rws;
    Rewriting_System * provisional_rws;
    General_Multiplier * gm;
    General_Multiplier * dgm;
    General_Multiplier * gm2;
    General_Multiplier * dgm2;
    FSA_Simple * wa;
    FSA_Simple * cosets;
    FSA_Simple * cclass;
    FSA_Simple * conjugator;
    FSA_Simple * diff2;
    FSA_Simple * minred;
    FSA_Simple * minkb;
    FSA_Simple * maxkb;
    FSA_Simple * rr;
    FSA_Simple * diff1c;
    FSA_Simple * diff2c;
    FSA_Simple * pdiff1;
    FSA_Simple * pdiff2;
    FSA_Simple * subwa;
    FSA_Buffer();
    APIMETHOD ~FSA_Buffer();
    APIMETHOD FSA * load(Container * container,String filename,unsigned ga_flags,
                         bool coset_system,MAF * maf);
};

/*****************************************************************************/

/* Class Presentation_Data is a home for various properties to be available
   read-only in public members. The author intensely dislikes getter/setter
   type methods. */
class Presentation_Data
{
  public:
    Presentation_Type presentation_type;
    Ordinal nr_generators; // for a coset system this includes the coset
                           // symbol and named subgroup generators and
                           // inverses
    Ordinal coset_symbol;  // this is INVALID_VALUE unless this is a coset
                           // system in which case it is also the number of
                           // G generators.
    Owned_String filename; // both this and next member are 0 unless the
    Owned_String name;     // Presentation is created from a GAP record
    Owned_String original_filename; /* for a coset system, the name of the group
                                       presentation */
    Owned_String subgroup_filename; /* for a coset system, the name of the
                                       subgroup file used to generate it */
    Total_Length max_relator_length; /* lhs+rhs for longest relator
                                        coset system relators are ignored */
    bool g_level;
    bool h_level;
    bool has_cancellation;
    bool inversion_difficult;/* true if we can expect to to be difficult
                                to recognise inverse pairs of reduced words */
    bool is_confluent;       /* this is always false unless we read in a
                                confluent rewriting system. This no longer
                                gets set to true if Rewriter_Machine succeeds
                                in building a confluent rewriting machine
                                from axioms */
    bool is_coset_system;
    bool is_normal_coset_system;
    /* The next two properties are set true if and when they are discovered.
       This is inconsistent with how is_confluent is now being handled */
    bool is_coset_finite; /* set true if the language of accepted words of
                             interest is finite. For a coset system this
                             means the language of words beginning with _H,
                             otherwise it is the entire accepted language */
    bool is_g_finite;   /* set true if the language of the entire rewriting
                           system is finite. This is a stronger property
                           than is_coset_finite. */
    bool is_group;/* only set true if inverses are specified for all generators.
                     So this being false, does not mean the object does NOT
                     describe a group */
    bool g_is_group; /* set true if inverses are specified for all the g generators */
    bool is_rubik;
    bool is_short; // true if reduction cannot increase word length
    bool is_shortlex; // true if order is true shortlex
  public:
    Presentation_Data(Presentation_Type presentation_type);
    virtual ~Presentation_Data() {};
    /*
       right_multipliers()
       sets start to the first valid right multiplier for a word ending with
       rvalue, and end to the last valid multiplier + 1 (so the returned
       valued are suitable for use in a C/C++ style for loop:
          for (g = start;g < end;g++)
       left_multipliers() returns a similar range for the symbols which
       can prefix a word beginning with lvalue
    */
    void right_multipliers(Ordinal * start,Ordinal * end,Ordinal rvalue) const
    {
      if (!is_coset_system || rvalue < 0)
      {
        *start = 0;
        *end = nr_generators;
      }
      else if (rvalue <= coset_symbol)
      {
        *start = 0;
        *end = coset_symbol;
      }
      else
      {
        *start = coset_symbol;
        *end = nr_generators;
      }
    }
    void left_multipliers(Ordinal * start,Ordinal * end,Ordinal lvalue) const
    {
      if (!is_coset_system)
      {
        *start = 0;
        *end = nr_generators;
      }
      else if (lvalue < coset_symbol)
      {
        *start = 0;
        *end = coset_symbol+1;
      }
      else
      {
        *start = coset_symbol+1;
        *end = nr_generators;
      }
    }
};

// flags for add_axiom
const unsigned AA_ADD_TO_RM = 1;
const unsigned AA_WARN_IF_REDUNDANT = 2;
const unsigned AA_DEDUCE = 4; // Make simple deductions:
                              // i.e. word-differences, conjugations etc.
const unsigned AA_ELIMINATE = 8;//attempt to prove this axiom is unnecessary
                                //before adding it by doing a more elaborate
                                //expansion of the current set of equations
const unsigned AA_POLISH = 16;  // instruction to Rewriter_Machine to add
                                // polished form of this axiom to the
                                // presentation
const unsigned AA_DEFAULT = AA_ADD_TO_RM|AA_WARN_IF_REDUNDANT|AA_POLISH;
const unsigned AA_RESUME = 32; // "axiom" is an equation deduced in an earlier run
/* flag AA_DEDUCE_INVERSE permits MAF to call set_inverse() if the equation
specified in add_axiom() is of the form x*y=IdWord and either x or y is not
currently invertible. Note that this means there is an assumption that
y*x=IdWord, which does not follow. So you should only pass this flag if you are
certain that the object you are dealing with is a group. MAF uses this flag
internally when determining a presentation for a subgroup using Schreier
generators. It might also be useful if a MAF object from something like
a .kbprog file, but one in which the inverses[] field is missing */
const unsigned AA_DEDUCE_INVERSE = 64;
const unsigned AA_INTERNAL = 128;

/* Class Presentation contains all the basic information given in a
   group or monoid presentation, or in a KBMAG style RWS file (but not
   the various KBMAG options that KBMAG records in such files), and just
   those methods which do not require any automata to have been grown.

   The author does not remember why he made this a separate class from MAF.
*/
class Presentation : protected Presentation_Data
{
  friend class Group_Automata;
  friend class RWS_Reader;
  friend class Subalgebra_Descriptor;
  protected:
    Ordinal *inverses;
    /* member group_words is only allocated when analysing a coset system
       with generator names.
       The nth entry is a word for the nth generator using only the main
       generators. So for the G generators it is just the word consisting
       of the generator alone. For the coset symbol _H it is the
       empty/padding word. For the H generators is the word given in the
       subalgebra descriptor, or the inverse of the group word for its
       inverse. */
    Packed_Word *group_words;
    Packed_Equation_List * axioms;
    Packed_Equation_List * polished_axioms;
    Generator_Permutation * permutation;
    bool container_ours;
    bool flags_set;
    Alphabet & real_alphabet;  // this must come before alphabet. We don't want
                               // to expose this except as a const
    const Alphabet * base_alphabet;
  public:
    const Alphabet & alphabet;
    Container & container;
    Presentation(Container & container,Alphabet_Type alphabet_type,
                 Presentation_Type presentation_type,bool container_ours);
    ~Presentation();
    APIMETHOD bool add_axiom(String lhs,String rhs,unsigned flags = AA_DEFAULT);
    APIMETHOD bool add_axiom(String lhs,String rhs,unsigned flags, Parse_Error_Handler &handler);
    APIMETHOD bool add_axiom(const Word &lhs,const Word &rhs,unsigned flags = AA_DEFAULT);
    APIMETHOD void add_permutation(String str);
    void delete_permutations();
    bool has_permutations() const
    {
      return permutation != 0;
    }
    APIMETHOD void set_generators(String gens);
    APIMETHOD void set_inverse(Ordinal g,Ordinal ig,unsigned flags = 0);
    APIMETHOD void set_inverses(String inverse);
    /* set_nr_generators()/set_next_generator() can be used where it is
       desired to add the generators one at a time.
       set_nr_generators() completely removes any information about
       generators already present, and must be called first.
       Individual generators can only be set consecutively */
    APIMETHOD bool set_next_generator(String glyph);
    APIMETHOD bool rename_generator(Ordinal g,String glyph);
    APIMETHOD bool set_nr_generators(Element_Count nr_letters_);
    // const methods from now on
    unsigned axiom_count(bool polished) const;
    unsigned axiom_count() const
    {
      return Presentation::axiom_count(polished_axioms!=0);
    }
    const Linked_Packed_Equation * first_axiom(bool polished) const;
    const Linked_Packed_Equation * first_axiom() const
    {
      return first_axiom(polished_axioms!=0);
    }
    Ordinal generator_count() const { return nr_generators;};
    APIMETHOD bool relators(Word_Collection * wc) const;
    // return a list of subgroup generators
    // NB this is taken from the coset equations in the cos file,
    // so might not match the sub file exactly.
    APIMETHOD bool subgroup_generators(Word_Collection * wc) const;
    const Alphabet & group_alphabet() const
    {
      /* Usually this is just the same as alphabet. But in a coset system
         it is the alphabet of the original group presentation, so less the
         coset symbol and subgroup generators. */
      return *base_alphabet;
    }
    APIMETHOD const Packed_Word & group_word(Ordinal g) const;
    void group_word(Ordinal_Word * answer, const Word & word) const;

    /* abelianise() returns the word that is the shortlex freeabelian
       reduction of the input word */
    APIMETHOD void abelianise(Word * rword,const Word & word) const;
    /* abelian_multiply returns the word that is the shortlex freeabelian
       reduction of the product of the words */
    APIMETHOD void abelian_multiply(Word *answer,const Word &w0,const Word & w1) const;
    /* abelian_invert() returns the word that is the shortlex freeabelian
       reduction of the inverse of the input word */
    APIMETHOD bool abelian_invert(Word * inverse_word,const Word & word) const;
    Ordinal inverse(Ordinal generator) const
    {
      return inverses[generator];
    }
    /* Neither of these invert methods return the reduced inverse. Instead
       they return the free group inverse, i.e. the word formed by reading
       from the right and writing from the left the inverse of each
       generator in the word */
    bool invert(Letter * iword,size_t *buffer_size,String start_word) const;
    APIMETHOD bool invert(Word * inverse_word,const Word & word,
                          Word_Length length = WHOLE_WORD) const;
    /* free_reduce() performs juset the reductions possible based on
       the inverse relations */
    APIMETHOD void free_reduce(Word * answer,const Word & word) const;

    // parse a GAP syntax expression into the specified word
    APIMETHOD bool parse(Ordinal_Word *word,String string,String_Length length = WHOLE_STRING) const;
    // parse a GAP syntax expression into a newly allocated word.
    // It is caller's responsibility to delete the word afterwards.
    APIMETHOD Ordinal_Word * parse(String string,
                                   String_Length length = WHOLE_STRING) const;
    APIMETHOD Ordinal_Word * parse(String string,
                                   String_Length length,
                                   Parse_Error_Handler &error_handler) const;
    const Presentation_Data & properties() const
    {
      return *this;
    }
    // Remaining methods are for internal use only
  public:
    void add_polished_axiom(Node_Manager &,const Word & lhs,const Word & rhs);
  private:
    virtual bool insert_axiom(const Word & lhs,const Word & rhs,
                              unsigned flags) = 0;
  protected:
    void set_flags();
};

/* WR_PREFIX_ONLY removes the last letter from the word before reduction
   and then puts it back afterwards.
*/
const unsigned WR_PREFIX_ONLY = 1;

/* If WR_ONCE is specified, reduce() will only improve the word once, so that
   the answer is not necessarily fully reduced. There are several reasons
   you might want to do this:

     It allows multiple RWS equations to be created from one Diff_Reduce
     reduction.

     If several Diff_Reduce buffers are available using WR_ONCE on the
     first and then using a second to fully reduce the word will improve
     performance when a tree of words is being searched.
*/
const unsigned WR_ONCE = 2;
/* WR_CHECK_ONLY checks for reducibility but not does perform the reduction */
const unsigned WR_CHECK_ONLY = 4;
/* WR_ASSUME_L2 tells the reducer to assume that no prefix of the word to
   be reduced is reducible. This may improve performance in some circumstances
*/
const unsigned WR_ASSUME_L2 = 8;
/* WR_SCHREIER_LABEL only applies to reduction of words using a coset MIDFA
   multiplier. If specified then the returned word is prefixed by the
   ordinals corresponding to the initial states used during the reduction,
   followed by the coset symbol, followed by the coset representative.
   So the returned word looks like the RHS of a coset equation in a coset
   system. WR_H_LABEL is similar, but uses the actual h word appearing in
   the label. */
const unsigned WR_SCHREIER_LABEL = 16;
const unsigned WR_H_LABEL = 32;

class Word_Reducer
{
  public:
    virtual ~Word_Reducer() {};
    /* The return value from reduce() is the number of separate reduction
       steps that were required to reduce the word. If WR_ONCE or
       WR_CHECK_ONLY are passed implementations can return 1 or 0 and
       are not expected to count the number of reductions that would be made
       if the word were to be fully reduced. */
    virtual unsigned reduce(Word * word,const Word & start_word,
                            unsigned flags = 0,const FSA * wa = 0) = 0;
    virtual bool reducible(const Word & word,Word_Length length = WHOLE_WORD);
    // this method is used internally in programs that construct automata when
    // reduction is through a provisional Word_Reducer. In such cases it may
    // be necessary for a Word_Reducer to be deleted and recreated.
    virtual bool broken() const
    {
      return false;
    }
};

/* logging flags */
const unsigned LOG_WORD_DIFFERENCES = 1;
const unsigned LOG_EQUATIONS = 2;
const unsigned LOG_DERIVATIONS = 4;
/* LOG_TRACE is set while a proof trace is being written (the -trace
   option). Like LOG_DERIVATIONS it makes MAF record how each equation is
   derived, but nothing is printed */
const unsigned LOG_TRACE = 8;
const unsigned LOG_RECORD = LOG_DERIVATIONS|LOG_TRACE;

/* flags for create_from_input(),create_rws() etc. */
const unsigned CFI_NAMED_H_GENERATORS = 1;
const unsigned CFI_NO_INVERSES = 2;
const unsigned CFI_CREATE_RM = 4;
const unsigned CFI_ALLOW_CREATE = 8;
const unsigned CFI_DEFAULT = CFI_NAMED_H_GENERATORS;
const unsigned CFI_CS_EXPECTED = 16; /* The input should be a coset system */
const unsigned CFI_RESUME = 32;
const unsigned CFI_RAW = 64;
const unsigned CFI_REQUIRE_NAMED_H_GENERATORS = 128;

/* -strategy flags in options.strategy */
const unsigned MSF_RESPECT_RIGHT = 1;
const unsigned MSF_LEFT_EXPAND = 2;
const unsigned MSF_AGGRESSIVE_DISCARD = 4;
const unsigned MSF_SELECTIVE_PROBE = 8;
const unsigned MSF_CONJUGATE_FIRST = 16;
const unsigned MSF_NO_FAVOUR_SHORT = 32;
const unsigned MSF_DEEP_RIGHT = 64;
const unsigned MSF_USE_ERAS = 128;
const unsigned MSF_CLOSED = 256;

const Byte BYTE_OPTION_UNSET = 100;

class MAF : public Presentation
{
  public:
    enum MAF_Strategy
    {
      MS_Naive, /* tries to make MAF behave like a very simple KB implementation */
      MS_Vanilla, /* tries to make MAF behave like a basic KB implementation */
      MS_Easy,    /* suitable only for "easy" presentations - for which it is fast */
      MS_Quick,   /* includes slightly more extra stuff, but still performs well */
      MS_Short,   /* usually best for difficult presentations in shortlex input files */
      MS_Balanced, /* very similar to MS_Short, but expands more equations sooner */
      MS_Wreath, /* usually best for difficult presentations using wreath orderings */
      MS_Sparse, /* strategy for automatic structures */
      MS_Era,   /* possible strategy for presentations like SL(2,p) balanced presentation */
      MS_Long /* probably best strategy for shortlex SL(2,p) balanced presentations */
    };
    struct Options
    {
      // log_flags enable additional diagnostics according to LOG_ flags above
      unsigned log_flags;
      unsigned log_level;
      // used by mafauto.cpp
      unsigned max_multiplier_attempts; // number of times to try building
                                       // multiplier after word_acceptor stabilises
      // options used by maf_rm.cpp and associated modules
      Element_Count max_equations;
      unsigned min_time;
      unsigned timeout;
      unsigned max_time;
      unsigned strategy;
      Byte filters;
      Byte probe_style;
      Byte expansion_order;
      Byte repeat;
      Byte balancing_flags;  // balancing options are used by maf_we.cpp to
      Byte left_balancing;   // control how equations are balanced
      Byte right_balancing;  //
      Byte secondary;
      Byte partial_reductions;
      Byte special_overlaps;
      Byte conjugation;
      Byte collapse;
      Byte detect_finite_index;
      bool assume_confluent;
      bool dense_rm; // true if sparse nodes are not to be used
      bool no_equations;
      // used by maf_jm.cpp
      bool no_early_repair;
      bool no_throttle;
      // used by maf_ew.cpp
      bool swap_bad;
      bool unbalance;
      bool ignore_h_length;
      bool no_h;
      // maf_rm.cpp only
      bool check_inverses;
      bool consider_secondary;
      bool weed_secondary;
      bool expand_all;
      bool extended_consider;
      bool no_deductions;
      bool no_pool;
      bool no_prune;
      bool fast;
      bool slow;
      bool tight;
      bool no_kb;
      bool force_differences;
      bool no_differences;
      bool differences;
      // used by maf_dt.cpp (mostly)
      bool no_half_differences;
      // used by mafauto.cpp
      bool emulate_kbmag;
      bool force_multiplier;
      bool is_kbprog;
      bool no_weak_acceptor;
      bool use_schreier_generators;
      bool validate;  // true if axiom check wanted
      bool validate_inverses; // true if axiom check should check inverse relations
      bool write_success;
      // used by subpres.cpp
      bool eliminate;
      bool no_composite;
      Total_Length pool_above;
      Total_Length max_overlap_length;
      Word_Length no_pool_below;
      Word_Length max_stored_length[2];
      Options() :
        log_flags(0),
        log_level(1),
        max_equations(0),
        max_multiplier_attempts(10),
        max_time(0),
        min_time(0),
        timeout(20),
        strategy(~0u),
        filters(7),
        probe_style(BYTE_OPTION_UNSET),
        expansion_order(0),
        secondary(BYTE_OPTION_UNSET),
        balancing_flags(3),
        left_balancing(BYTE_OPTION_UNSET),
        right_balancing(BYTE_OPTION_UNSET),
        collapse(0),
        conjugation(1),
        partial_reductions(1),
        repeat(0),
        special_overlaps(1),
        assume_confluent(false),
        check_inverses(false),
        consider_secondary(false),
        dense_rm(false),
        detect_finite_index(false),
        eliminate(false),
        emulate_kbmag(false),
        expand_all(false),
        extended_consider(false),
        fast(false),
        force_differences(false),
        force_multiplier(false),
        ignore_h_length(false),
        is_kbprog(false),
        no_composite(true),
        no_differences(false),
        differences(false),
        no_deductions(false),
        no_early_repair(false),
        no_equations(false),
        no_h(false),
        no_half_differences(false),
        no_kb(false),
        no_pool(false),
        no_prune(false),
        no_throttle(false),
        no_weak_acceptor(false),
        unbalance(true),
        slow(false),
        swap_bad(false),
        tight(false),
        use_schreier_generators(false),
        validate(false),
        validate_inverses(false),
        weed_secondary(false),
        write_success(false),
        pool_above(0),
        no_pool_below(0),
        max_overlap_length(0)
      {
        max_stored_length[0] = max_stored_length[1] = 0;
      }
      void set_strategy(MAF_Strategy predefined_strategy)
      {
        switch (predefined_strategy)
        {
          case MS_Naive:
            unbalance = false;
            no_pool = true;
            filters = 0;
            partial_reductions = 0;
            special_overlaps = 0;
            conjugation = 0;
            secondary = 0;
            right_balancing = 0;
            left_balancing = 0;
            strategy = MSF_RESPECT_RIGHT|MSF_DEEP_RIGHT;
            no_half_differences = no_deductions = true;
            probe_style = 0;
            expansion_order = 4; //to show how bad it can get!
            break;
          case MS_Vanilla:
            unbalance = false;
            no_pool = true;
            filters = 0;
            partial_reductions = 0;
            special_overlaps = 0;
            conjugation = 0;
            secondary = 0;
            strategy = MSF_RESPECT_RIGHT|MSF_USE_ERAS;
            no_half_differences = no_deductions = true;
            probe_style = 0;
            expansion_order = 2; // 4 would be more "naive" but renders this strategy useless
            break;
          case MS_Easy:
            partial_reductions = 0;
            special_overlaps = 0;
            conjugation = 0;
            secondary = 0;
            no_half_differences = no_deductions = true;
            strategy = MSF_RESPECT_RIGHT|MSF_AGGRESSIVE_DISCARD;
            probe_style = 0;
            expansion_order = 1;
            break;
          case MS_Quick:
            partial_reductions = 1;
            special_overlaps = 0;
            conjugation = 1;
            secondary = 1;
            strategy = MSF_RESPECT_RIGHT|MSF_AGGRESSIVE_DISCARD;
            probe_style = 0;
            expansion_order = 1;
            break;
          case MS_Short:
            partial_reductions = 1;
            special_overlaps = 1;
            conjugation = 1;
            secondary = 1;
            strategy = MSF_RESPECT_RIGHT|MSF_AGGRESSIVE_DISCARD;
            probe_style = 3;
            expansion_order = 1;
            break;
          case MS_Balanced:
            partial_reductions = 1;
            special_overlaps = 1;
            conjugation = 1;
            secondary = 1;
            strategy = MSF_RESPECT_RIGHT|MSF_AGGRESSIVE_DISCARD|MSF_LEFT_EXPAND;
            probe_style = 3;
            break;
          case MS_Wreath:
            partial_reductions = 1;
            special_overlaps = 1;
            conjugation = 1;
            secondary = 3;
            strategy = MSF_RESPECT_RIGHT;
            probe_style = 3;
            expansion_order = 2;
            break;
          case MS_Sparse:
            partial_reductions = 1;
            special_overlaps = 1;
            conjugation = 1;
            secondary = 2;
            strategy = MSF_RESPECT_RIGHT|MSF_AGGRESSIVE_DISCARD|MSF_SELECTIVE_PROBE;
            probe_style = 2;
            expansion_order = 1;
            break;
          case MS_Long:
            partial_reductions = 1;
            special_overlaps = 0;
            conjugation = 0;
            secondary = 1;
            strategy = MSF_RESPECT_RIGHT|MSF_NO_FAVOUR_SHORT|MSF_CLOSED;
            probe_style = 0;
            expansion_order = 4;
            break;
          case MS_Era:
            expand_all = true;
            partial_reductions = 1;
            probe_style = 3;
            special_overlaps = 0;
            conjugation = 0;
            expansion_order = 4;
            strategy = MSF_RESPECT_RIGHT|MSF_NO_FAVOUR_SHORT|MSF_CLOSED|MSF_USE_ERAS;
            break;
        }
      }
    };
  private:
    FSA_Buffer real_fsas;
    FSA_Buffer real_group_fsas;
    Rewriter_Machine *rm;
    Diff_Reduce * validator;
    Word_Reducer * wr;
    mutable Word_Reducer * provisional_wr;
    /* FSAs */
    Group_Automata * automata;
  public:
    bool aborting;
    Options options;
    const FSA_Buffer &fsas;
  public:
    static Container * create_container(Platform * platform = 0);
    // create an empty MAF object
    static MAF * create(Container * container = 0,
                        const Options * options = 0,
                        Alphabet_Type = AT_String,
                        Presentation_Type = PT_General);
    // options to Create the MAF object from KBMAG records
    // create_from_rws() can create coset systems or ordinary MAF systems
    // the file contents determine what is created
    static MAF * create_from_rws(String filename,Container * container = 0,
                                 unsigned flags = CFI_DEFAULT,
                                 const Options * options = 0);
    // create_from_substructure() creates a coset system from an existing rws
    // and substructure file, or else opens a previously created coset system
    static MAF * create_from_substructure(String group_filename,
                                          String subgroup_suffix = 0,
                                          Container * container = 0,
                                          unsigned flags =
                                          CFI_DEFAULT|CFI_CS_EXPECTED,
                                          const Options * options = 0);
    // create_from_input() can create a MAF object either from file group_name
    // (which might be a coset system!) or or it can
    // use create_from_substructure() to first create and then open a
    // newly created coset system.
    static MAF * create_from_input(bool want_coset_system,
                                   String group_filename,
                                   String subgroup_suffix = 0,
                                   Container * container = 0,
                                   unsigned flags = CFI_DEFAULT,
                                   const Options * options = 0)
    {
      return want_coset_system ?
        create_from_substructure(group_filename,subgroup_suffix,container,flags,options) :
        create_from_rws(group_filename,container,flags & CFI_CREATE_RM+CFI_RESUME,options);
    }
    /* deprecated method used by author's package Spirofractal */
    static FSA_Simple * translate_acceptor(const MAF & maf_new,
                                           const MAF & maf_start,
                                           const FSA & fsa_start,
                                           const Word_List &wl);

    /* translate_acceptor translates the word acceptor for a group into
       a new alphabet. A word acceptor for the group must already have
       been computed. The object which computes the translation must be
       a coset system with named subgenerators and inverses and the
       subgroup must have index 1, and the coset general multiplier must
       have been computed */
    APIMETHOD FSA_Simple * translate_acceptor();
    ~MAF();

    void delete_fsa(Group_Automaton_Type ga_type);

    void give_up()
    {
      aborting = true;
    }

    // Todd Coxeter methods. Enumeration is always against the relators of
    // the underlying group. For the first method subgroup generators must
    // be specified explicitly. For the second method MAF will create the
    // appropriate set of subgroup generators based on the type of MAF object
    // If you pass a non-zero value for coset_table then this is returned to
    // you and you own the object and must delete it.
    APIMETHOD Language_Size enumerate_cosets(FSA **coset_table,
                                             const Word_Collection &normal_closure_generators,
                                             const Word_Collection &subgroup_generators,
                                             const TC_Enumeration_Options & options);
    APIMETHOD Language_Size enumerate_cosets(FSA **coset_table,const TC_Enumeration_Options & options);

    /* realise_rm() can be used to convert a MAF instance created without
       a Rewriter_Machine instance into one with one */
    APIMETHOD void realise_rm();

    /* grow_automata() works as follows:
       1) As many of the automata indicated by the value
          (save_flags | retain_flags)
          that can be built are built. Each of these flags is first anded
          with ~exclude_flags if it has the value GA_ALL.
       2) The automata that exist and are specified in save_flags are written
          out to text files using suitable names and suffixes.
       3) If a non-null FSA_Buffer * is specified any automata indicated by
          retain_flags are assigned to the appropriate part of the buffer.
          The caller now owns these and is responsible for deleting them.
          If you prefer not to have to do this then specify 0 for retain_flags
          and reload the FSAs you require afterwards with load_fsas() or
          load_group_fsas().
       4) All remaining automata are deleted apart from any that MAF retains
          for its own purposes.
       A non-zero value for exclude_flags should only be specified if at
       least one of the other flags is GA_ALL, otherwise it will be ignored.
    */
    APIMETHOD void grow_automata(FSA_Buffer * buffer = 0,
                                 unsigned save_flags = GA_ALL,
                                 unsigned retain_flags = 0,
                                 unsigned exclude_flags = 0);
    APIMETHOD void grow_geodesic_automata(bool need_multiplier,
                                          unsigned fsa_format_flags);
    APIMETHOD void import_difference_machine(const FSA_Simple & dm);
    const FSA * load_fsas(unsigned flags)
    {
      // The caller does not own the FSAs returned by this function
      // they will be deleted when the MAF object is deleted
      // the return value is the last FSA loaded
      return real_fsas.load(&container,filename,flags,is_coset_system,this);
    }
    const FSA * load_group_fsas(unsigned flags)
    {
      // The caller does not own the FSAs returned by this function
      // they will be deleted when the MAF object is deleted
      // the return value is the last FSA loaded
      if (!is_coset_system)
        return load_fsas(flags);
      return real_group_fsas.load(&container,original_filename,flags,false,this);
    }
    APIMETHOD bool load_reduction_method(Group_Automaton_Type flag);

    Word_Reducer * take_word_reducer()
    {
      /* method to allow the caller to take over management of a Word_Reducer
         from MAF. The caller is now responsible for deleting the object
         (if it exists!) */
      Word_Reducer * answer = wr;
      wr = 0;
      return answer;
    }
    APIMETHOD bool output_gap_presentation(Output_Stream * stream,
                                           bool change_alphabet) const;
    APIMETHOD void print(Output_Stream * stream) const;
    APIMETHOD void save(String filename);
    APIMETHOD void save_as(String basefilename,String suffix="")
    {
      String_Buffer sb;
      save(String::make_filename(&sb,"",basefilename,suffix));
    }
#ifdef FSA_INCLUDED
    void save_fsa(FSA * fsa,String suffix,unsigned format_flags = 0) const
    {
      fsa->save_as(filename,name,suffix,format_flags);
    }
#endif
    APIMETHOD void save_fsa(FSA * fsa,Group_Automaton_Type gat,unsigned format_flags = 0) const;
    APIMETHOD void set_validation_fsa(FSA_Simple * fsa);

    const FSA_Buffer &group_fsas() const
    {
      if (!is_coset_system)
        return real_fsas;
      return real_group_fsas;
    }
    // is_valid_equation is intended() for internal validation only
    // it only works if set_validation_fsa() was called previously, or
    // if the is_rubik flag is set.
    // the return value is false if the elements corresponding to the words
    // are detected as being unequal, in which case error is assigned to
    // appropriate error text.
    APIMETHOD bool is_valid_equation(String *error,const Word & lhs,const Word & rhs) const;

    APIMETHOD String multiplier_name(String_Buffer *sb,const Word & word,bool coset_multiplier) const;
    bool reduction_available() const
    {
      return wr != 0;
    }

    Rewriter_Machine & rewriter_machine()
    {
      realise_rm();
      return *rm;
    }

    // method to parse a list of words contained in a GAP list
    APIMETHOD void read_word_list(Word_List *wl,String filename) const;
    APIMETHOD bool reduce(String_Buffer * rword,String word) const;
    APIMETHOD unsigned reduce(Word * rword,const Word & word,unsigned flags = 0) const;
    APIMETHOD FSA_Simple *translate_acceptor(const MAF & maf_start,
                                             const FSA & fsa_start,
                                             const Word_List &xlat_wl) const;
    Alphabet * get_alphabet()
    {
      return rm ? 0 : &real_alphabet;
    }

    // methods for processing words

    /* Check a word for reducibility */
    bool reducible(String word,String_Length length = WHOLE_STRING) const;
    bool reducible(const Word &word,Word_Length length = WHOLE_WORD) const;
    bool L2_accepted(String word,String_Length length = WHOLE_STRING) const;
    bool L2_accepted(const Word &word) const;
    bool polish_equation(Simple_Equation * se) const;
    APIMETHOD unsigned long order(Word *test,const FSA * wa) const;

    // methods for building automata
    /* mafauto.cpp */
    FSA_Simple * polish_difference_machine(FSA *old_diff,
                                           bool delete_old);
    FSA_Simple * tidy_difference_machine(FSA_Simple * fsa_) const;
    APIMETHOD FSA_Simple * build_acceptor_from_dm(const FSA *dm,
                                                  bool create_equations,
                                                  bool geodesic = false,
                                                  bool force_group_acceptor = false,
                                                  Word_Length max_context = 0,
                                                  State_Count max_states = 0);
    APIMETHOD FSA_Simple * build_acceptor_from_coset_table(const FSA & coset_table);
    APIMETHOD FSA_Simple * labelled_product_fsa(const FSA &fsa,
                                                bool accept_only = false) const;
    /* maf_mult.cpp */
    APIMETHOD bool check_axioms(const General_Multiplier &gm,
                                Compositor_Algorithm algorithm = CA_Hybrid,
                                bool check_inverses = false) const;
    /* mafcoset.cpp */
    APIMETHOD FSA_Simple * coset_table(Word_Reducer & coset_wr,
                                        const FSA & wa) const;
    void label_coset_table(FSA_Common * coset_table) const;
    void subgroup_generators_from_coset_table(Word_Collection * wc,
                                                const FSA &coset_table) const;
    /* mafconj.cpp */
    APIMETHOD FSA_Simple * build_class_table(const Rewriting_System &rws) const;
    APIMETHOD FSA_Simple * build_conjugator(const FSA_Simple &class_table,
                                           const Rewriting_System &rws) const;

    /* subpres.cpp */
    APIMETHOD MAF *subgroup_presentation(General_Multiplier & coset_gm,
                                         bool use_schreier_generators);
    APIMETHOD MAF *subgroup_presentation(const Rewriting_System & rws,
                                         bool use_schreier_generators) const;
    APIMETHOD MAF *rs_presentation(const FSA & coset_table) const;
    /* maf_subwa.cpp */
    APIMETHOD FSA_Simple * subgroup_word_acceptor(Word_Reducer &coset_wr,
                                              const FSA * coset_dm = 0);

  private:
    MAF(Container & container,
        Alphabet_Type alphabet_type,Presentation_Type presentation_type,
        bool container_ours,const Options * options);
    // virtual method from Presentation we have to implement
    bool insert_axiom(const Word & lhs,const Word & rhs,unsigned flags);

    FSA_Simple * reducer(const General_Multiplier & multiplier,
                        const FSA & difference_machine,
                        bool full_length = false);

    void subgroup_relators(const General_Multiplier &gm,
                           Word_List * relator_wl,
                           bool use_schreier_generators);
};

#ifdef FSA_INCLUDED
// Multiplier deals with all FSAs which accept pairs of accepted words (u,v)
// and in which the accepting state has a label w that make uw==v a true
// equation
class Multiplier : public Delegated_FSA
{
  private:
    struct Node;
    struct Node_List;
    Label_ID * containing_label;
    Sorted_Word_List & multipliers;
    Element_Count nr_multipliers;
    mutable Word_Length max_depth;
    mutable Word_Length valid_length;
    mutable Node_List * stack;
  public:
    Multiplier(const Multiplier & other);
    Multiplier(FSA &other,bool owner_);
    ~Multiplier();
    /* composite() builds a multiplier supporting the list of multipliers
       specified in new_multipliers. The method should only be called on a
       General_Multiplier instance, or a Multiplier whose multipliers are
       a superset of the generators, or to extract the multiplier for a
       subset (usually of one element) of the multipliers of a previously
       built multiplier. The method may fail, even when in theory it should
       not, if asked to create a composite multiplier for some subgroup of a
       subgroup, because it may not recognise which composites it requires
       (because of reductions). In future I may add a parameter to allow
       a composite to be merged in with a gm. In this case success
       would be guaranteed, and performance would be improved in the case
       where you want to build such a multiplier.
    */

    APIMETHOD FSA_Simple * composite(const MAF & maf,
                                     const Word_Collection &new_multipliers) const;
    FSA_Simple * composite(const MAF & maf,
                           const Word & word) const;
    APIMETHOD FSA_Simple * kbmag_multiplier(const MAF & maf,
                                  const Word & word,
                                  const String prefix = "_x");
    bool labels_consistent() const
    {
      /* if containing_label is non-zero then each multiplier
         appears in exactly one of the labels associated with
         an accept state. This means that if the multiplier accepts
         (u,v) for m1 and (u,v) for m2 and (x,y) for m1 it also
         accepts (x,y) for m2. If not then it is possible that
         (x,y) might not be accepted for m2.
         This will happen only when the multiplier is a
         determinised coset multiplier.
         This information might be useful when axiom checking,
         and is useful when doing reduction since if labels are
         consitent then any words except the empty word in the label
         for the initial state are relators (e.g. trivial genrators) */
      return containing_label!=0;
    }
    bool is_multiplier(const Word & word,Element_ID * multiplier_nr = 0) const;
    bool label_contains(Label_ID label,Element_ID multiplier_nr) const;
    bool label_contains(Label_ID label,const Word & word) const
    {
      Element_ID multiplier_nr;
      if (is_multiplier(word,&multiplier_nr))
        return label_contains(label,multiplier_nr);
      return false;
    }
    /* On entry start_word must be accepted. It is multiplied by the
       multiplier contained in the specified position of multipliers
       and the answer is put in answer. The return code is true if
       the multiplication was successful and false if not, which can
       happen if any of the following are true:
       1) start_word is not accepted by the acceptor corresponding to the multiplier
       2) multiplier_nr is not one of the valid multipliers
       3) the multiplier FSA is invalid
    */
    Element_Count multiplier_count() const
    {
      return nr_multipliers;
    }
    const Word * multiplier(Element_ID word_nr) const;
    const Sorted_Word_List &multiplier_list() const
    {
      return multipliers;
    }
    bool multiply(Ordinal_Word * answer,const Word & start_word,
                  Element_ID multiplier_nr,State_ID * initial_state = 0) const;
    // Compute the multiplier for the inverses of the elements in this
    // multiplier
    FSA_Simple *inverse(const MAF &maf) const;
  private:
    void set_multipliers();
};

// General_Multiplier is a class for the specific Multiplier that
// has the identity and all the generators as its multipliers
// This multiplier can also do reduction
class General_Multiplier : public Multiplier
{
  private:
    Element_ID * multiplier_nr;
  public:
    General_Multiplier(const General_Multiplier & other);
    General_Multiplier(FSA_Simple &other,bool owner_);
    ~General_Multiplier();
    bool label_for_generator(Label_ID label,Ordinal g) const;
    unsigned reduce(Word * answer,const Word & start_word,unsigned flags = 0) const;
  private:
    void set_multiplier_nrs();
};

class GM_Reducer : public Word_Reducer
{
  public:
    General_Multiplier & gm;
    GM_Reducer(General_Multiplier & gm_) :
      gm(gm_)
    {}
    unsigned reduce(Word * word,const Word & start_word,unsigned flags,
                    const FSA *)
    {
      return gm.reduce(word,start_word,flags);
    }
};
#endif

#endif

