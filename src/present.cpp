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


// $Log: present.cpp $
// Revision 1.16  2010/06/10 13:57:46Z  Alun
// All tabs removed again
// Revision 1.15  2010/06/08 06:33:54Z  Alun
// Flag setting handling changed to make it more robust.
// relators() and rename_generator() method added
// Revision 1.14  2009/10/13 23:30:26Z  Alun
// Can no longer rely on input errors bombing out.
// Revision 1.13  2009/09/12 18:47:57Z  Alun
// printf() style functions cleaned up considerably
// code now compiles cleanly in g++ both 32-bit and 64-bit
// all tabs removed from source files
// Revision 1.12  2008/12/27 18:14:30Z  Alun
// Support for normal subgroups, abelianisation and free reduction added
// Revision 1.12  2008/11/02 18:59:24Z  Alun
// Early November checkin
// Completely reworked pool and Difference_Tracker
// Ported to Darwin
// Revision 1.11  2008/10/09 11:25:02Z  Alun
// Need to call set_flags() in case where there is no Rewriter_Machine.
// Moved checking of coset system axioms
// Revision 1.10  2008/09/26 07:24:50Z  Alun
// Final version built using Indexer.
// Revision 1.5  2007/12/20 23:25:42Z  Alun
//

#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include "maf.h"
#include "mafword.h"
#include "equation.h"
#include "container.h"
#include "maf_trace.h"

Presentation_Data::Presentation_Data(Presentation_Type presentation_type_) :
  nr_generators(0),
  presentation_type(presentation_type_),
  is_confluent(false),
  is_g_finite(false),
  is_coset_finite(false),
  is_group(false),
  g_is_group(false),
  is_rubik(presentation_type_ == PT_Rubik),
  is_coset_system(false),
  is_normal_coset_system(false),
  is_short(false),
  is_shortlex(false),
  has_cancellation(false),
  max_relator_length(0),
  inversion_difficult(true)
{
}

/**/

Presentation::Presentation(Container & container_,Alphabet_Type alphabet_type,
                           Presentation_Type Presentation_Type_,bool container_ours_) :
  Presentation_Data(Presentation_Type_),
  container(container_),
  container_ours(container_ours_),
  permutation(0),
  real_alphabet(*Alphabet::create(alphabet_type,container_)),
  alphabet(real_alphabet),
  base_alphabet(&real_alphabet),
  inverses(0),
  group_words(0),
  axioms(new Packed_Equation_List),
  polished_axioms(0),
  flags_set(false)
{
  real_alphabet.attach();
  base_alphabet->attach();
}

/**/

Presentation::~Presentation()
{
  if (inverses)
    delete [] inverses;
  if (group_words)
    delete [] group_words;
  base_alphabet->detach();
  delete axioms;
  if (polished_axioms)
    delete polished_axioms;
  delete_permutations();
  real_alphabet.detach();
  if (container_ours)
    delete &container;
}

/**/

unsigned Presentation::axiom_count(bool polished) const
{
  return polished && polished_axioms ? polished_axioms->count() : axioms->count();
}

const Linked_Packed_Equation * Presentation::first_axiom(bool polished) const
{
  return polished && polished_axioms ? polished_axioms->first() : axioms->first();
}

/**/

bool Presentation::subgroup_generators(Word_Collection * wc) const
{
  Ordinal_Word word(group_alphabet());
  wc->empty();
  if (!is_coset_system)
    return false;
  for (const Linked_Packed_Equation * axiom = first_axiom();
       axiom;
       axiom = axiom->get_next())
  {
    Simple_Equation se(alphabet,*axiom);
    if (se.lhs_word.value(0) == coset_symbol)
    {
      word = Subword(se.lhs_word,1);
      Word_Length r = se.rhs_word.length();
      Ordinal g;
      while (r-- > 0 && (g = se.rhs_word.value(r)) < coset_symbol)
        word.append(inverse(g));
      wc->add(word);
    }
  }
  return true;
}

/**/

bool Presentation::relators(Word_Collection * wc) const
{
  Ordinal_Word word(group_alphabet());
  wc->empty();
  bool retcode = true; // BUG FIXME this code clearly looks wrong; retcode is not returned
  const Ordinal nr_generators = group_alphabet().letter_count();

  for (const Linked_Packed_Equation * axiom = first_axiom();
       axiom;
       axiom = axiom->get_next())
  {
    Simple_Equation se(alphabet,*axiom);
    if (se.last_letter() <  nr_generators) {
      if (se.relator(&word,*this))
        wc->add(word);
      else
        retcode = false;
    }
  }
  return true;
}

/**/

bool Presentation::parse(Ordinal_Word *word,String string,String_Length length) const
{
  return alphabet.parse(word,string,length,this,container);
}

/**/

Ordinal_Word * Presentation::parse(String string,String_Length length,Parse_Error_Handler &error_handler) const
{
  return alphabet.parse(string,length,this,error_handler);
}

/**/

Ordinal_Word * Presentation::parse(String string,String_Length length) const
{
  return parse(string,length,container);
}

/**/

void Presentation::add_permutation(String str)
{
  Generator_Permutation * p = new Generator_Permutation;
  /* Permutations can be used to reduce the number of axioms in
     groups with symmetric generators */
  p->next = permutation;
  p->word = parse(str);
  if (p->word)
    permutation = p;
  else
    delete p;
}

/**/

void Presentation::delete_permutations()
{
  while (permutation)
  {
    Generator_Permutation * p = permutation->next;
    delete permutation;
    permutation = p;
  }
}

/**/

void Presentation::add_polished_axiom(Node_Manager &,
                                      const Word & lhs_word,
                                      const Word & rhs_word)
{
  if (!polished_axioms)
    polished_axioms = new(Packed_Equation_List);
  polished_axioms->add(lhs_word,rhs_word);
}

/**/

bool Presentation::add_axiom(const Word &lhs_word,const Word &rhs_word,unsigned flags)
{
  bool retcode = true;
  Word_Length l = lhs_word.length();
  Word_Length r = rhs_word.length();

  if (!l)
  {
    /* ensure lhs word is not empty */
    if (!r)
      return false;
    return add_axiom(rhs_word,lhs_word,flags);
  }

  if (flags & AA_DEDUCE_INVERSE && r == 0)
  {
    if (l==2)
    {
      Ordinal x = lhs_word.value(0);
      Ordinal X = lhs_word.value(1);
      if (inverse(x) == INVALID_SYMBOL || inverse(X) == INVALID_SYMBOL)
      {
        set_inverse(x,X,AA_INTERNAL);
        set_flags();
      }
    }
    else if (l==1)
    {
      /* In this case the generator is trivial, but we need to make sure
         it has an inverse otherwise simplify won't work */
      Ordinal x = lhs_word.value(0);
      if (inverse(x) == INVALID_SYMBOL)
      {
        set_inverse(x,x,AA_INTERNAL);
        set_flags();
      }
    }
  }

  if (!flags_set)
    set_flags();

  /* These tests used to be performed after the axiom was inserted, which
     is clearly undesirable. However, because they have moved here
     set_flags() has to be careful not to destroy group_words if set_flags
     has been called before.
  */
  if (is_coset_system)
  {
    /* We check that any equation being inserted into a coset system obeys
       certain rules.

       First we divide words into four types.
       1. The empty word.
       2. G-words, containing only G-generators,
       3. H-words, containing only H-generators,
       4. mixed words, words containing both types of generator and words
          containing the coset symbol.

       A mixed word is only valid if it has the form u * coset_symbol * v
       where u is either empty or an H-word and v is either empty or a G-word.

       Now we can state the rules.

       1) Both sides must be valid words
       2) If one side of the equation is a G-word or an H-word then the other
       side must be a word of the same type or the empty word.
       3) If one side of an equation is a mixed word, the other side must
       also be a mixed word. These equations are the "coset equations".

       Formerly the code also required that the LHS of an equation
       involving mixed words had an empty u part. This is true of all the
       coset equations created when a coset system is generated, but might
       not be true in an output file if either H-generators were named
       but not all invertible, or the underlying object is a monoid and
       not a group. So this check has been removed.
    */

    String_Buffer sb1;
    String_Buffer sb2;
    String lhs = lhs_word.format(&sb1);
    String rhs = rhs_word.format(&sb2);
    const Ordinal * values = lhs_word.buffer();
    Word_Length i;
    Ordinal lvalue = values[0];


    if (lvalue < coset_symbol)
    {
      /* if the LHS word starts with a G-generator both sides must be
         G-words, except that the RHS might be empty */
      for (i = 1; i < l ; i++)
        if (values[i] >= coset_symbol)
        {
          container.input_error("LHS has incorrect form for a coset system."
                                " LHS beginning with %s should not\ncontain "
                                "generator %s. Offending equation is:\n"
                                "%s=%s\n",alphabet.glyph(lvalue).string(),
                                alphabet.glyph(values[i]).string(),lhs.string(),rhs.string());
          return false;
        }
      values = rhs_word.buffer();
      for (i = 0; i < r ; i++)
        if (values[i] >= coset_symbol)
        {
          container.input_error("RHS has incorrect form for a coset system."
                                " If LHS begins with %s RHS should not\n"
                                "contain generator %s. Offending equation"
                                " is:\n%s=%s\n",alphabet.glyph(lvalue).string(),
                                alphabet.glyph(values[i]).string(),lhs.string(),
                                rhs.string());
          return false;
        }
    }
    else if (values[l-1] > coset_symbol)
    {
      /* If the LHS ends with an H-generator it must be an H-word.
         Check both sides are H-words (but allowing an empty RHS) */
      lvalue = values[l-1];
      for (i = 0; i < l ; i++)
        if (values[i] <= coset_symbol)
        {
          container.input_error("LHS has incorrect form for a coset system."
                                " LHS ending with %s should not\ncontain "
                                "generator %s. Offending equation is:\n"
                                "%s=%s\n",alphabet.glyph(lvalue).string(),
                                alphabet.glyph(values[i]).string(),
                                lhs.string(),rhs.string());
          return false;
        }
      values = rhs_word.buffer();
      for (i = 0; i < r ; i++)
        if (values[i] <= coset_symbol)
        {
          container.input_error("RHS has incorrect form for a coset system."
                                " If LHS end with %s RHS should not\n"
                                "contain generator %s. Offending equation"
                                " is:\n%s=%s\n",alphabet.glyph(lvalue).string(),
                                alphabet.glyph(values[i]).string(),
                                lhs.string(),rhs.string());
          return false;
        }
    }
    else /* The LHS began with _H or an h generator, but ended with _H or a
            G-generator, so both sides must have the form u * coset_symbol*v */
    {
      bool coset_matched = false;
      for (i = 0; i < l ; i++)
        if (values[i] == coset_symbol && !coset_matched)
          coset_matched = true;
        else if (values[i] < coset_symbol && !coset_matched)
        {
          container.input_error("LHS has incorrect form for a coset system."
                                " If LHS begins with %s it should not\n"
                                "contain generator %s before %s. Offending "
                                "equation is:\n%s=%s\n",
                                alphabet.glyph(lvalue).string(),
                                alphabet.glyph(values[i]).string(),
                                alphabet.glyph(coset_symbol).string(),
                                lhs.string(),rhs.string());
          return false;
        }
        else if (values[i] >= coset_symbol && coset_matched)
        {
          container.input_error("LHS has incorrect form for a coset system."
                                " If LHS begins with %s it should not\n"
                                "contain generator %s after %s. Offending "
                                "equation is:\n%s=%s\n",
                                alphabet.glyph(lvalue).string(),
                                alphabet.glyph(values[i]).string(),
                                alphabet.glyph(coset_symbol).string(),
                                lhs.string(),rhs.string());
          return false;
        }

      // coset_matched is now true otherwise we would have returned an error

      values = rhs_word.buffer();
      coset_matched = false;
      for (i = 0; i < r ; i++)
        if (values[i] == coset_symbol && !coset_matched)
          coset_matched = true;
        else if (values[i] < coset_symbol && !coset_matched)
        {
          container.input_error("RHS has incorrect form for a coset system."
                                " If LHS begins with %s RHS should not\n"
                                "contain generator %s before %s. Offending "
                                "equation is:\n%s=%s\n",
                                alphabet.glyph(lvalue).string(),
                                alphabet.glyph(values[i]).string(),
                                alphabet.glyph(coset_symbol).string(),
                                lhs.string(),rhs.string());
          return false;
        }
        else if (values[i] >= coset_symbol && coset_matched)
        {
          container.input_error("RHS has incorrect form for a coset system."
                                " If LHS begins with %s RHS should not\n"
                                "contain generator %s after %s. Offending "
                                "equation is:\n%s=%s\n",
                                alphabet.glyph(lvalue).string(),
                                alphabet.glyph(values[i]).string(),
                                alphabet.glyph(coset_symbol).string(),
                                lhs.string(),rhs.string());
          return false;
        }

      if (!coset_matched)
      {
        container.input_error("RHS has incorrect form for a coset system."
                              " If LHS begins with %s RHS should \n"
                              "contain coset symbol %s. Offending "
                              "equation is:\n%s=%s\n",
                              alphabet.glyph(lvalue).string(),
                              alphabet.glyph(coset_symbol).string(),
                              lhs.string(),rhs.string());
        return false;
      }

      if (presentation_type >= PT_Coset_System_With_Generators && r == 2 &&
          values[1] == coset_symbol)
      {
        Subword gword((Word &) lhs_word,1,WHOLE_WORD);
        Ordinal hvalue = values[0];
        Ordinal_Word iword(alphabet);
        group_words[hvalue] = gword;
        if (inverse(hvalue) != INVALID_SYMBOL)
        {
          invert(&iword,gword,l-1);
          group_words[inverse(hvalue)] = iword;
        }
      }
    }
  }

  if (proof_trace)
    proof_trace->given_axiom(lhs_word,rhs_word);
  if (flags & AA_ADD_TO_RM)
    retcode = insert_axiom(lhs_word,rhs_word,flags);
  else
  {
    axioms->add(lhs_word,rhs_word);
  }

  if (lhs_word.value(0) < group_alphabet().letter_count())
  {
    Total_Length tl = lhs_word.length()+rhs_word.length();
    if (tl > max_relator_length)
      max_relator_length = tl;
  }

  return retcode;
}

/**/

bool Presentation::add_axiom(String lhs,String rhs,unsigned flags,
                             Parse_Error_Handler &error_handler)
{
  Ordinal_Word * lhs_word = parse(lhs,WHOLE_STRING,error_handler);
  Ordinal_Word * rhs_word = parse(rhs,WHOLE_STRING,error_handler);
  bool retcode = add_axiom(*lhs_word,*rhs_word,flags);
  delete lhs_word;
  delete rhs_word;
  return retcode;
}

/**/

bool Presentation::add_axiom(String lhs,String rhs,unsigned flags)
{
  return add_axiom(lhs,rhs,flags,container);
}

/**/

void Presentation::set_generators(String identifiers_)
{
  real_alphabet.set_letters(identifiers_);
  nr_generators = alphabet.letter_count();
  if (inverses)
    delete inverses;
  inverses = new Ordinal[nr_generators+1];
  Ordinal j;
  for (j = 0; j < nr_generators;j++)
    inverses[j] = INVALID_SYMBOL;
  inverses[j] = nr_generators;
  flags_set = false;
}

/**/

bool Presentation::set_nr_generators(Element_Count nr_generators_)
{
  bool retcode = real_alphabet.set_nr_letters(nr_generators = Ordinal(nr_generators_));
  if (inverses)
    delete inverses;
  inverses = new Ordinal[nr_generators+1];
  Ordinal j;
  for (j = 0; j < nr_generators;j++)
    inverses[j] = INVALID_SYMBOL;
  inverses[j] = nr_generators;
  flags_set = false;
  return retcode;
}

/**/

bool Presentation::set_next_generator(Glyph glyph)
{
  return real_alphabet.set_next_letter(glyph);
}

/**/

bool Presentation::rename_generator(Ordinal g,Glyph glyph)
{
  return real_alphabet.rename_letter(g,glyph);
}

/**/

void Presentation::set_inverses(String inverses_)
{
  Ordinal_Word s(alphabet,2);
  Ordinal_Word id_word(alphabet,(Word_Length)0);
  Ordinal_Word & inverse_word = *parse(inverses_);
  const Ordinal * values = inverse_word.buffer();
  Word_Length l = inverse_word.length();
  Word_Length i;

  if (l & 1)
  {
    container.input_error("Inverse string passed to add_inverses must be even in length\n");
    return;
  }

  for (i = 0;i < l;i += 2)
    set_inverse(values[i],values[i+1],AA_INTERNAL);
  flags_set = false;
  delete &inverse_word;
}

/**/

void Presentation::set_flags()
{
  /* This is a critical method, but the user cannot call it directly.
     It usually only needs to be called once, just before the first
     axiom is inserted. However, if AA_DEDUCE_INVERSE is being used it
     may need to be called repeatedly.
     Formerly this method used to be called "just in case" a lot, but
     this could cause the group_words member to go wrong.
     So now the code does not recreate this member if set_flags() ought
     to have set it correctly before.
  */

  bool inverse_found = false;

  coset_symbol = INVALID_SYMBOL;
  Ordinal nr_non_invertible = 0;
  Ordinal first_non_invertible = INVALID_SYMBOL;
  Ordinal last_non_invertible = INVALID_SYMBOL;
  bool set_normal = false;
  g_is_group = is_group = true;
  for (Ordinal j = 0;j < nr_generators;j++)
  {
    String s = alphabet.glyph(j);
    if (strcmp(s,"_H")==0)
      coset_symbol = j;
    else if (strcmp(s,"_N")==0)
    {
      coset_symbol = j;
      set_normal = true;
    }

    if (inverses[j] == INVALID_SYMBOL)
    {
      is_group = false;
      nr_non_invertible++;
      if (first_non_invertible == INVALID_SYMBOL)
        first_non_invertible = j;
      last_non_invertible = j;
    }
    else
      inverse_found = true;
  }

  if (presentation_type == PT_General && nr_non_invertible > 0 &&
      (alphabet.order_type() == WO_Wreath_Product ||
       alphabet.order_type() == WO_Right_Wreath_Product ||
       alphabet.order_type() == WO_Coset) &&
      coset_symbol != INVALID_SYMBOL)
  {
    if (coset_symbol == nr_generators -1)
    {
      presentation_type = PT_Simple_Coset_System;
      container.progress(1,"Presentation being treated as simple coset system\n");
    }
    else if (last_non_invertible == coset_symbol)
    {
      presentation_type = PT_Coset_System_With_Inverses;
      container.progress(1,"Presentation being treated as coset system with"
                           " sub generators and inverses\n");
    }
    else
    {
      presentation_type = PT_Coset_System_With_Generators;
      container.progress(1,"Presentation being treated as coset system with "
                           "sub generators\n");
    }
  }
  if (presentation_type >= PT_Simple_Coset_System)
  {
    bool valid = false;
    for (;;) /* Coded as a loop to allow us to use break as a goto end! */
    {
      if (alphabet.order_type() != WO_Wreath_Product &&
          alphabet.order_type() != WO_Right_Wreath_Product &&
          alphabet.order_type() != WO_Coset)
      {
        container.error_output("Coset systems must use a Wreath Product type ordering\n");
        break;
      }
      if (nr_non_invertible == 0)
      {
        container.error_output("There is no symbol that can act as the coset"
                               " separator\n");
        break;
      }

      if (coset_symbol == INVALID_SYMBOL)
      {
        coset_symbol = first_non_invertible;
        container.error_output("The special _H symbol was not found. ");
        if (coset_symbol)
          container.error_output("Generator %s is being assumed to be the"
                                 " coset symbol instead.\n",
                                 alphabet.glyph(coset_symbol).string());
        else
        {
          container.error_output("\n");
          break;
        }
      }

      if (coset_symbol == 0)
      {
        container.error_output("There are no main object generators!\n");
        break;
      }

      if (presentation_type > PT_Simple_Coset_System &&
          coset_symbol == nr_generators - 1)
      {
        presentation_type = PT_Simple_Coset_System;
        container.error_output("No named subgroup generatators are specified"
                               " after the coset symbol, so the\nsystem will "
                               "be treated as a simple coset system instead.\n");
      }
      if (presentation_type == PT_Coset_System_With_Inverses &&
          last_non_invertible > coset_symbol)
      {
        presentation_type = PT_Coset_System_With_Generators;
        container.error_output("One or more of the subgroup generators has no"
                               " inverse so the system will be\n treated as a"
                               " coset system with generators instead\n");
      }
      valid = true;
      break;
    }

    if (!valid)
    {
      container.error_output("The rules for coset systems have not been "
                             "followed so the the presentation\nwill be "
                             "treated generically.\n");
      presentation_type = PT_General;
    }
    else
    {
      is_coset_system = true;
      if (presentation_type >= PT_Coset_System_With_Generators)
      {
        if (!flags_set)
        {
          if (group_words)
            delete [] group_words;
          group_words = new Packed_Word[nr_generators];
          Ordinal_Word ow(alphabet,1);
          Ordinal g;
          for (g = 0;g < coset_symbol;g++)
          {
            ow.set_code(0,g);
            group_words[g] = ow;
          }
          ow.set_length(0);
          for (;g < nr_generators;g++)
            group_words[g] = ow;
        }
      }
      else
        is_normal_coset_system = set_normal;
      real_alphabet.set_coset_order(coset_symbol);
      base_alphabet->detach();
      base_alphabet = alphabet.truncated_alphabet(coset_symbol);
      base_alphabet->attach();
      if (!original_filename || !subgroup_filename)
      {
        /* In fact we always get here, because maf_sub.cpp can't set these
           fields until after the coset system object has been created.
           But it may be that the user specifies an explicit coset system
           filename.
           We try to work out the name of the corresponding group and
           substructure files.
           We don't report any problems here, because it is possible that
           we have been called from maf_sub.cpp, and that it knows the
           right names to use*/
        if (filename != 0)
        {
          String suffix = filename.suffix();
          if (*suffix)
          {
            size_t l = suffix.length();
            if (l >= 4)
            {
              int sub = 0;
              if (String(suffix + l-4).is_equal("_cos"))
                sub |= 2;
              else if (suffix.matching_length(".cos") == 4)
                sub |= 1;
              if (sub == 3)
                sub = 0; /* The coset suffix looks like .cosxxx_cos so
                            we can't tell whether the substructure was
                            .subxxx_cos or .cosxxx. The user deserves
                            to suffer.  */
              original_filename = filename.clone_substring(0,suffix - filename);
              if (sub != 0)
              {
                Letter * temp_filename = sub==1 ?
                                         filename.clone() :
                                         filename.clone_substring(0,filename.length()-4);
                if (sub == 1)
                {
                  size_t i = String(temp_filename).suffix() - temp_filename;
                  memcpy(temp_filename+i,".sub",4);
                }
                subgroup_filename = temp_filename;
              }
            }
          }
        }
      }
    }
  }
  if (is_coset_system ? first_non_invertible < coset_symbol : !is_group)
    g_is_group = false;

  if (!is_group && !is_coset_system)
  {
    if (presentation_type == PT_Group || presentation_type == PT_Rubik)
      if (inverse_found)
        container.progress(1,"Not all generators have generator inverses, so type"
                           " has been changed to\nmonoid\n");
    presentation_type = PT_General;
  }
  if (inverse_found && presentation_type == PT_Monoid_Without_Cancellation)
    presentation_type = PT_General;
  if (!inverse_found)
  {
    if (presentation_type == PT_Group || presentation_type == PT_Rubik)
      container.progress(1,"No generators are known to be invertible so type"
                         " has been changed to monoid\nwithout cancellation\n");
    presentation_type = PT_Monoid_Without_Cancellation;
  }

  if (is_group)
  {
    has_cancellation = true;
    if (presentation_type != PT_Rubik)
      presentation_type = PT_Group;
  }
  if (presentation_type == PT_Monoid_With_Cancellation)
    has_cancellation = true;
  if (coset_symbol == INVALID_SYMBOL)
    coset_symbol = nr_generators;
  is_shortlex = alphabet.order_is_shortlex();
  is_short = is_shortlex ||
             alphabet.order_is_geodesic() ||
             alphabet.order_is_effectively_shortlex();
  /* we are setting inversion_difficult for all orders that can
     increase word length indefinitely. But it would probably be worth
     making this optional. Difficulty increases with the perversity of the
     ordering, with orderings that are not ordinally similar to 1,2,3,...
     definitely perverse, and wreath type orders with multiple levels
     particularly so, since even when each generator has the same level as
     its inverse, the maximal subword will typically move to the opposite
     side from where the orderings wants it, so that the free inverse is
     likely to be reducible and require a lot of reductions, which will
     most likely increase the length alarmingly. */
  inversion_difficult = !group_alphabet().order_is_geodesic() &&
                        !group_alphabet().order_is_effectively_shortlex() &&
                        group_alphabet().order_needs_moderation();
  if (is_coset_system)
  {
    g_level = h_level = true;
    Ordinal g = 0;
    unsigned level = alphabet.level(g);
    for (g = 1; g < coset_symbol;g++)
      if (alphabet.level(g) != level)
        g_level = false;

    level = alphabet.level(g);
    for (; g < nr_generators;g++)
      if (alphabet.level(g) != level)
        h_level = false;
  }
  else
  {
    g_level = is_short;
    h_level = true;
  }
  flags_set = true;
}

/**/

void Presentation::set_inverse(Ordinal g,Ordinal ig,unsigned flags)
{
  Ordinal_Word s(alphabet,2);
  Ordinal_Word id_word(alphabet);
  bool needed = false;

  if (inverses[g] == INVALID_SYMBOL || ig < inverses[g])
  {
    inverses[g] = ig;
    needed = true;
  }
  if (ig != INVALID_SYMBOL && (inverses[ig] == INVALID_SYMBOL || g < inverses[ig]))
  {
    inverses[ig] = g;
    needed = true;
  }

  if (needed)
  {
    if (max_relator_length < 2)
      max_relator_length = 2;
    s.set_code(0,g);
    s.set_code(1,ig);
    if (flags & AA_POLISH)
      axioms->add(s,id_word);
    if (g != ig)
    {
      s.set_code(0,ig);
      s.set_code(1,g);
      if (flags & AA_POLISH)
        axioms->add(s,id_word);
    }
  }
}

/**/

bool Presentation::invert(Word * inverse_word,const Word & word,Word_Length length) const
{
  if (length==WHOLE_WORD)
    length = word.length();
  inverse_word->allocate(length,true); // in case we are inverting in place
  const Ordinal * old_ordinals = word.buffer();
  Ordinal * new_ordinals = inverse_word->buffer();
  for (Word_Length i = 0; i < length--;i++)
  {
    Ordinal g1 = inverse(old_ordinals[i]);
    Ordinal g2 = inverse(old_ordinals[length]);
    if (g1 == INVALID_SYMBOL || g2 == INVALID_SYMBOL)
      return false;
    new_ordinals[i] = g2;
    new_ordinals[length] = g1;
  }
  return true;
}

/**/

bool Presentation::invert(Letter * iword,size_t * buf_size,String start_word) const
{
  bool retcode = false;
  Ordinal_Word *parsed = parse(start_word);
  if (parsed)
  {
    Ordinal_Word answer(alphabet,parsed->length());
    if (invert(&answer,*parsed,parsed->length()))
    {
      *buf_size = answer.format(iword,*buf_size);
      retcode = true;
    }
    else
      *buf_size = 0;
    delete parsed;
  }
  return retcode;
}

const Packed_Word & Presentation::group_word(Ordinal g) const
{
  return group_words[g];
}

/**/

void Presentation::group_word(Ordinal_Word * answer,const Word & start_word) const
{
  /* Converts start_word into a word in the main generators by removing the
     coset separator symbol and all subgroup generators.
  */
  Word_Length i = 0;

  if (answer == &start_word)
  {
    Ordinal_Word clone(start_word);
    group_word(answer,clone);
    return;
  }

  answer->set_length(0);
  if (is_coset_system)
  {
    const Ordinal * values = start_word.buffer();
    for (i = 0; values[i] > coset_symbol;i++)
      *answer += Ordinal_Word(*base_alphabet,group_words[values[i]]);
    if (values[i] == coset_symbol)
      i++;
  }
  *answer += Subword((Word &) start_word,i,WHOLE_WORD);
}

/**/

void Presentation::abelianise(Word * answer,const Word & word) const
{
  /* All the methods to do with abelianisation arguably belong in MAF.
     But since we don't need access to any of its methods or data it
     seems sensible to put them here. In any case, these methods are
     mainly used when manipulating presentations, */

  Word_Length l = word.length();
  if (answer->buffer() != word.buffer())
    *answer = word;
  Ordinal * values = answer->buffer();

  /* See if the word is already in its correct form */
  bool ok = true;
  for (Word_Length i = 0; i+1 < l;i++)
    if (values[i+1] < values[i])
    {
      ok = false;
      break;
    }
    else
    {
      Ordinal ig = inverse(values[i]);
      if (ig != INVALID_SYMBOL)
      {
        if (ig+1 < values[i] || values[i] < ig+1 ||
            (ig >= values[i] && values[i+1]==ig))
        {
          ok = false;
          break;
        }
      }
    }

  if (ok)
    return;
  Word_Length position = 0;
  while (position < l)
  {
    Ordinal g = values[position];
    Word_Length i = position;
    while (++i < l)
    {
      if (values[i] == inverse(g))
      {
        l -= 2;
        if (position < l)
        {
          if (i != l)
          {
            g = values[position] = values[l];
            values[i] = values[l+1];
          }
          else
            g = values[position] = values[l+1];
          i = position;
        }
      }
    }
    position++;
  }
  answer->set_length(l);
  position = 0;
  while (position < l)
  {
    Ordinal best_g = values[position];
    Word_Length best_position = position;
    Word_Length i = position;
    while (++i < l)
    {
      if (values[i] < best_g)
      {
        best_g = values[i];
        best_position = i;
      }
    }
    values[best_position] = values[position];
    values[position] = best_g;
    position++;

    while (position < l && values[position] == best_g)
      position++;
  }
}

/**/

void Presentation::abelian_multiply(Word * answer,const Word &w0, const Word & w1) const
{
  if (answer->buffer() == w0.buffer())
  {
    Ordinal_Word temp(w0);
    abelian_multiply(answer,temp,w1);
    return;
  }
  if (answer->buffer() == w1.buffer())
  {
    Ordinal_Word temp(w1);
    abelian_multiply(answer,w0,temp);
    return;
  }
  Word_Length l0 = w0.length();
  Word_Length l1 = w1.length();
  Word_Length p0 = 0;
  Word_Length p1 = 0;
  const Ordinal * w0_values = w0.buffer();
  const Ordinal * w1_values = w1.buffer();

  answer->allocate(l0+l1,false);
  Ordinal *to = answer->buffer();

  while (p0 < l0 && p1 < l1)
  {
    if (w0_values[p0] == inverse(w1_values[p1]))
    {
      p0++;
      p1++;
    }
    else if (w0_values[p0] == w1_values[p1])
    {
      *to++ = w0_values[p0++];
      *to++ = w1_values[p1++];
    }
    else if (w0_values[p0] < w1_values[p1])
      *to++ = w0_values[p0++];
    else
      *to++ = w1_values[p1++];
  }
  while (p0 < l0)
    *to++ = w0_values[p0++];
  while (p1 < l1)
    *to++ = w1_values[p1++];
  answer->set_length(to-answer->buffer());

  abelianise(answer,*answer);
}

/**/

bool Presentation::abelian_invert(Word * answer,const Word & word) const
{
  Word_Length l = word.length();
  if (answer->buffer() != word.buffer())
    *answer = word;
  Ordinal * values = answer->buffer();
  for (Word_Length i = 0;i < l;i++)
  {
    Ordinal ig = inverse(values[i]);
    if (ig == INVALID_SYMBOL)
      return false;
    values[i] = ig;
  }
  abelianise(answer,*answer);
  return true;
}

/**/

void Presentation::free_reduce(Word * answer,const Word &w0) const
{
  Word_Length old_length = w0.length();
  Word_Length new_length = 0;
  Ordinal last_g = PADDING_SYMBOL;

  if (answer->buffer() != w0.buffer())
    *answer = w0;

  Ordinal * values = answer->buffer();

  for (Word_Length i = 0; i < old_length;i++)
  {
    if (inverse(values[i]) != last_g)
      values[new_length++] = last_g = values[i];
    else
    {
      new_length--;
      last_g = new_length ? values[new_length-1] : PADDING_SYMBOL;
    }
  }
  answer->set_length(new_length);
}
