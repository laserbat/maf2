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


// $Log: maf_ew.cpp $
// Revision 1.2  2010/06/10 13:57:37Z  Alun
// All tabs removed again
// Revision 1.1  2010/04/03 07:40:30Z  Alun
// New file.
//

#include "maf_rm.h"
#include "mafnode.h"
#include "maf_nm.h"
#include "container.h"
#include "maf_wdb.h"
#include "maf_ew.h"

Equation_Word::Equation_Word(Node_Manager &nm_,
                             const Word & start_word) :
  nm(nm_),
  allocated(start_word.length()),
  language(language_A),
  valid_length(0),
  state(new Node_ID[allocated+1]),
  values_base(new Ordinal[allocated+1])
{
  allocate(start_word.length(),false);
  unpack(start_word);
}

/**/

Equation_Word::Equation_Word(Node_Manager &nm_,Node_Handle state_,Ordinal transition) :
  nm(nm_),
  language(state_->language(nm_)),
  valid_length(state_->valid_length(nm_)),
  allocated((word_length = state_->length(nm_)) + (transition>=0)),
  state(new Node_ID[allocated+1]),
  values_base(values = new Ordinal[allocated+1]),
  certificate(nm_.current_certificate())
{
  state[0] = nm.start();
  state_->read_states_and_values(values,state+1,nm,state_,word_length);
  if (transition >= 0)
    append(transition);
  else
    values[word_length] = TERMINATOR_SYMBOL;
#ifdef DEBUG
  format(&debug_word);
#endif
}

/**/

Equation_Word::Equation_Word(Node_Manager &nm_,Word_Length length_) :
  nm(nm_),
  allocated(length_),
  language(language_A),
  valid_length(0),
  state(new Node_ID[length_+1]),
  values_base(values = new Ordinal[length_+1])
{
  set_length(0);
}

/**/

Equation_Word::Equation_Word(const Equation_Word & other) :
  Word(other), // needed to shut up g++ - does nothing
  nm(other.nm),
  allocated(other.word_length),
  state(new Node_ID[allocated+1]),
  values_base(values = new Ordinal[allocated+1])
{
  operator=(other); // This is going to be Word::operator=
}

/**/

Equation_Word & Equation_Word::operator=(Equation_Word & other)
{
  allocate(other.word_length,&other == this);
  if (other.values != values_base)
  {
    memcpy(values_base,other.values,word_length*sizeof(Ordinal));
    values = values_base;
  }
  if (&other != this)
  {
    certificate = other.certificate;
    valid_length = other.valid_length;
    if (valid_length)
      memcpy(state,other.state,(valid_length+1)*sizeof(Node_ID));
    language = other.language;
    reductions.merge(&other.reductions);
  }
#ifdef DEBUG
  format(&debug_word);
#endif
  return *this;
}

/**/

Equation_Word::~Equation_Word()
{
  delete [] state;
  delete [] values_base;
}

/**/

void Equation_Word::allocate(Word_Length new_length,bool keep)
{
  state[0] = nm.start();
  if (new_length > MAX_WORD)
    MAF_INTERNAL_ERROR(nm.maf.container,("Attempt to create overlong word\n"));

  if (new_length > allocated)
  {
    int i = ((new_length + 31) & ~31);
    if (i > MAX_WORD)
      allocated = MAX_WORD;
    else
      allocated = Word_Length(i);
    Node_ID * new_state = new Node_ID[allocated+1];
    Ordinal * new_values = new Ordinal[allocated+1];
    if (keep)
    {
      memcpy(new_state,state,(1+valid_length)*sizeof(Node_ID));
      memcpy(new_values,values,word_length*sizeof(Ordinal));
    }
    else
    {
      valid_length = 0;
      language = language_A;
      new_state[0] = nm.start();
    }
    delete [] values_base;
    values = values_base = new_values;
    delete [] state;
    state = new_state;
  }
  if (values != values_base)
  {
    if (keep)
      memcpy(values_base,values,word_length*sizeof(Ordinal));
    values = values_base;
  }
#ifdef DEBUG
  bool reformat = keep && new_length < word_length;
#endif
  values[word_length = new_length] = TERMINATOR_SYMBOL;
  if (valid_length > new_length)
    valid_length = new_length;
#ifdef DEBUG
  if (reformat)
    format(&debug_word);
#endif
}

/**/

void Equation_Word::read_word(Node_Handle node,bool read_state)
{
  allocate(node->length(nm),false);
  values = values_base;
  certificate = nm.current_certificate();
  if (read_state)
  {
    valid_length = node->valid_length(nm);
    state[0] = nm.start();
    node->read_states_and_values(values,state+1,nm,node,word_length);
  }
  else
  {
    valid_length = 0;
    node->read_values(values,nm,word_length);
  }
#ifdef DEBUG
  format(&debug_word);
#endif
  language = node->language(nm);
}

/**/

bool Equation_Word::read_inverse(Node_Reference nr,bool read_state,bool remember)
{
  Node_Reference inverse(0,0);

  if (remember)
  {
    nr = nr->reduced_node(nm,nr);
    if (!nr->is_final())
      inverse = nr->node(nm).set_inverse(nm,nr);

    if (inverse)
    {
      read_word(inverse,read_state);
      return true;
    }
  }

  if (!nr->read_inverse(this,nm))
    return false;
  language = language_A;
  valid_length = 0;
  reduce(0);
  return true;
}

/**/

Node_Reference Equation_Word::left_half_difference(Node_Handle state,
                                                   Ordinal g,bool create)
{
  /* by calling this twice we can calculate Xwy where w is the word for state */

  Node_Reference nr = state->left_half_difference(nm,state,g);
  if (nr.is_null() || !nr->is_valid())
  {
    read_inverse(state,true,true);
    if (g != PADDING_SYMBOL)
    {
      append(g);
      reduce();
    }
    nr = state_word(true);
    if (nr.is_null() && create)
      nr = get_node(word_length,0,0);
  }
  return nr;
}

/**/

Node_Reference Equation_Word::right_half_difference(Node_Handle state,Ordinal g,bool create)
{
  /* by calling this twice we can calculate Xwy where w is the word for state.
     It is slower than left_half_difference() because we need to read the
     word twice. */

  Node_Reference nr = state->right_half_difference(nm,state,g);
  if (!nr.is_null() && nr->is_valid())
    return nr;
  read_word(state,true);
  if (g != PADDING_SYMBOL)
    append(g);
  reduce();
  nr = state_word(true);
  if (nr.is_null() && create)
    nr = get_node(word_length,0,0);
  if (!nr.is_null())
    return nr->node(nm).set_inverse(nm,nr,Node_Reference(0,0),true);
  return nr;
}

/**/

Node_Reference Equation_Word::difference(State state,Ordinal g1, Ordinal g2,bool create)
{
  Word_Length l = state->length(nm);
  if (g1 != PADDING_SYMBOL)
    l++;
  if (g2 != PADDING_SYMBOL)
    l++;
  allocate(l,false);
  Ordinal * v = values;
  if (g1 != PADDING_SYMBOL)
    *v++ = nm.inverse(g1);
  state->read_values(v,nm);
  if (g2 != PADDING_SYMBOL)
    values[l-1] = g2;
  valid_length = 0;
  bool failed = false;
  return create ? realise_state() : (reduce(0,&failed),state_word(true));
}

Node_Reference Equation_Word::realise_state()
{
  /* reduce the word and ensure it has a full length state */
  bool failed = false;
  reduce(0,&failed);
  if (!failed)
  {
    Node_Reference nr = state_word(true);
    if (nr.is_null())
      nr = get_node(word_length,0,0);
    return nr;
  }
  return Node_Reference(0,0);
}

/**/

Node_Reference Equation_Word::re_realise_state(Node_Handle state,bool attach)
{
  if (!state->is_final())
    return state;
  read_word(state,true);
  if (attach)
  {
    Node_Reference nr = realise_state();
    state->detach(nm,state);
    if (!nr.is_null())
      nr->attach(nm);
    return nr;
  }
  else
    return realise_state();
}

/**/

void Equation_Word::append(Ordinal rvalue)
{
  Node_Reference nr(0,0);
  validate();
  if (valid_length == word_length)
  {
    nr = Node_Reference(nm,state[word_length])->transition(nm,rvalue);
    if (!nr->is_final())
    {
      valid_length++;
      language = language_L0;
      certificate = nm.current_certificate();
    }
    else
      language = nr->language(nm);
  }
  else if (language == language_L0)
    language = language_L3;
  else
    language = language_A;
  set_length(word_length+1);
  values[word_length-1] = rvalue;
  state[word_length] = nr;
#ifdef DEBUG
  format(&debug_word);
#endif
}

/**/

void Equation_Word::pad(Word_Length new_length)
{
  Word_Length save_length = word_length;
  if (values + new_length > values_base+allocated+1)
  {
    set_length(new_length-1);
    values[word_length = save_length] = TERMINATOR_SYMBOL;
  }
  for (Word_Length i = save_length+1; i < new_length;i++)
    values[i] = PADDING_SYMBOL;
}

/**/

void Equation_Word::join(const Equation_Word & first_word,
                         const Equation_Word & second_word,Word_Length offset)
{
  bool keep = this == &first_word || this==&second_word;
  Word_Length first_length = offset == WHOLE_WORD ? first_word.word_length : offset;
  Word_Length second_length = second_word.word_length;
  if (first_length + second_length > MAX_WORD)
    MAF_INTERNAL_ERROR(nm.maf.container,("Maximum Word_Length exceeded in Equation_Word::join()\n"));
  allocate(first_length+second_length,keep);
  if (!keep)
  {
    memcpy(values,first_word.values,first_length*sizeof(Ordinal));
    memcpy(values+first_length,
           second_word.values,second_length*sizeof(Ordinal));
  }
  else
  {
    if (this == &second_word)
      memmove(values+first_length,second_word.values,
              second_length*sizeof(Ordinal));
    if (this != &first_word)
      memcpy(values,first_word.values,first_length*sizeof(Ordinal));
    if (this != &second_word)
      memcpy(values+first_length,
           second_word.values,second_length*sizeof(Ordinal));
  }
  if (&first_word != this)
    valid_length = 0;
  certificate = Certificate();
  language = language_A;
#ifdef DEBUG
  format(&debug_word);
#endif
}

/**/

Equation_Word & Equation_Word::operator=(const Word &word)
{
  allocate(word.length(),false);
  unpack(word);
  return *this;
}

/**/

void Equation_Word::unpack(const Word & word)
{
  values = values_base;
  valid_length = 0;
  language = language_A;
  word_copy(values,word.buffer(),word_length);
#ifdef DEBUG
  format(&debug_word);
#endif
}

/**/

bool Equation_Word::reduce(unsigned flags,bool * failed)
{
  /* Convert word to reduced form by applying equations.
     Return value is false if word is reduced, or true if it gets reduced
     If prefix_only is true then we only reduce the prefix of the word.
     We do this when the word is supposed to be the LHS of an equation */

  bool reduction = false;
  bool ok = true;
  Node_Reference s;
  Ordinal rvalue = PADDING_SYMBOL;

  /* What gets recorded is the reductions of this call, which the caller
     passes on to the equation's derivation if it keeps the result. Anything
     left from an earlier call was a trial reduction of words since replaced */
  if (nm.maf.options.log_flags & LOG_RECORD)
    reductions.empty();

  if (nm.maf.reduction_available())
  {
    /* Most reduction in the automata growing parts of maf eventually gets to
       Equation_Word::reduce(). But some programs need to create a
       Rewriter_Machine after the automata have already been created. In this
       case the Rewriter_Machine reduction is probably not very good, and
       we shall get the best reduction by using maf's currently loaded reduction
       mechanism. However we still need to reduce the word after this, to
       ensure the state stack is correct. */
    if (nm.maf.reduce((Word *) this,*this,flags & AE_KEEP_LHS+AE_SECONDARY ? WR_PREFIX_ONLY : 0))
      valid_length = 0;
  }

  state[0] = nm.start();

  if (!word_length)
    flags = 0;
  else if (flags & AE_KEEP_LHS+AE_SECONDARY)
  {
    word_length--;
    rvalue = values[word_length];
    if (valid_length > word_length)
      valid_length = word_length;
  }
  validate();
  s = Node_Reference(nm,state[valid_length]);
  if (nm.rm.presentation_trivial())
  {
    reduction = word_length != 0;  // we must set return code properly!
    set_length(0);
  }

  while (valid_length < word_length)
  {
    s = s->transition(nm,values[valid_length]);
    state[++valid_length] = s;
    if (s->is_equation())
    {
      Node_Reference rhs = s->fast_reduced_node(nm);
      Word_Length rhs_length = rhs->length(nm);
      Word_Length lhs_length = s->lhs_length(nm);
      reduction = true;
      if (rhs_length > lhs_length)
      {
        /* This can happen when we are not using shortlex order */
        if (nm.pd.is_short)
          MAF_INTERNAL_ERROR(nm.maf.container,("Bad reduction encountered!\n"));
        reduction = ok = hard_reduce(flags,failed);
        if (ok)
          s = Node_Reference(nm,state[word_length]);
        break;
      }
      valid_length -= lhs_length;
      if (nm.maf.options.log_flags & LOG_RECORD)
        reductions.add_state(s,valid_length);
      if (valid_length || rhs->is_final())
        rhs->read_values(values+valid_length,nm,rhs_length);
      else
        rhs->read_states_and_values(values,state+1,nm,rhs,rhs_length);
      if (rhs_length < lhs_length)
      {
        memcpy(values+valid_length+rhs_length,values+valid_length+lhs_length,
               (word_length-valid_length-lhs_length)*sizeof(Ordinal));
        word_length -= lhs_length - rhs_length;
        values[word_length] = TERMINATOR_SYMBOL;
      }
#ifdef DEBUG_REDUCE
      format(&debug_word);
#endif
      if (!valid_length && !rhs->is_final())
        valid_length += rhs_length;
      s = Node_Reference(nm,state[valid_length]);
    }
  }

  if (ok)
  {
    language = language_L0;
    valid_length = word_length;
    certificate = nm.current_certificate();
  }
  if (flags & AE_KEEP_LHS+AE_SECONDARY)
  {
    set_length(word_length+1);
    values[word_length-1] = rvalue;
    if (ok)
    {
      state[word_length] = s = s->transition(nm,rvalue);
      if (!s->is_final())
        valid_length = word_length;
      else
      {
        Word_Length lhs_length = s->lhs_length(nm);
        if (lhs_length < word_length)
          language = lhs_length == valid_length && s->fast_is_primary() ? language_L2 : language_L3;
        else if (s->fast_is_primary())
          language = language_L1;
        else
          language = s->fast_primary(nm)->lhs_length(nm) == valid_length ? language_L2 : language_L3;
      }
    }
  }
#ifdef DEBUG
  format(&debug_word);
#endif
  return reduction;
}

/**/

bool Equation_Word::hard_reduce(unsigned flags,bool * failed)
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
  Node_Reference s(nm,state[valid_length]);
  Node_Reference rhs = s->fast_reduced_node(nm);

  if (valid_length == word_length && s->lhs_length(nm) == word_length &&
      !rhs->is_final())
  {
    /* In this case we can do an immediate reduction and return */
    if (nm.maf.options.log_flags & LOG_RECORD)
      reductions.add_state(s,0);
    allocate(rhs->length(nm));
    rhs->read_states_and_values(values,state+1,nm,rhs,word_length);
    valid_length = word_length;
    return true;
  }

  bool overflow = false;
  bool reducing = true;

  Ordinal_Word save(nm.maf.alphabet);
  Ordinal_Word stack_word(nm.maf.alphabet,MAX_WORD);
  Ordinal * stack = stack_word.buffer();
  if (failed)
    save = *this;

  Word_Length max_length = MAX_WORD;
  Word_Length start_length = word_length;
  unsigned long nr_reductions = 0;
  bool reported = false;
  unsigned long time = 1;
  if (flags & AE_KEEP_LHS+AE_SECONDARY)
    max_length--;
  Word_Length stack_top = 0;
  for (Word_Length i = word_length; i > valid_length;)
    stack[stack_top++] = values[--i];
  Total_Length total_length = word_length;
  if (values != values_base)
    set_length(valid_length);
  for (;;)
  {
    if (reducing && s->is_final())
    {
      Node_Reference rhs = s->fast_reduced_node(nm);
      Word_Length rhs_length = rhs->length(nm);
      Word_Length lhs_length = s->lhs_length(nm);
      /* we deliberately don't immediately detect situation where maximum word
         length has been exceeded, because we may recover in the end. We
         don't let the length exceed the maximum for long however, as that
         probably is a waste of time */
      total_length -= lhs_length - rhs_length;
      if (lhs_length==valid_length && !rhs->is_final())
      {
        /* In this case we can do an immediate reduction */
        if (nm.maf.options.log_flags & LOG_RECORD)
          reductions.add_state(s,0);
        allocate(rhs_length);
        rhs->read_states_and_values(values,state+1,nm,rhs,rhs_length);
        valid_length = word_length;
      }
      else
      {
        valid_length -= lhs_length;
        if (nm.maf.options.log_flags & LOG_RECORD)
          reductions.add_state(s,valid_length);
        word_length = valid_length;
        if (stack_top + rhs_length > MAX_WORD)
        {
          overflow = true;
          break;
        }
        rhs->read_reversed_values(stack+stack_top,nm);
        stack_top += rhs_length;
      }
      s = Node_Reference(nm,state[valid_length]);
      nr_reductions++;
      if (!(char) nr_reductions)
      {
        if (total_length > max_length)
        {
          overflow = true;
          break;
        }

        if (nm.rm.critical_status(2,reported ? 1 : 3,"Performing a difficult"
                                  " reduction. Reductions %lu. Length %d/%d.\n",
                                  nr_reductions,word_length,total_length))

        {
          time++;
          reported = true;
          if (failed && time >= 10)
          {
            nm.rm.critical_status(2,0,"Aborting reduction\n");
            /* In this case we won't restore the original value
               because in a difficult collapse it may help to
               keep where we have got to so far */
            reducing = false;
            *failed = true;
          }
        }
      }
    }

    if (stack_top)
    {
      if (word_length >= allocated)
      {
        Word_Length save = word_length;
        if (total_length > max_length)
        {
          overflow = true;
          break;
        }
        allocate(Word_Length(total_length),true);
        word_length = save;
      }
      Ordinal g = values[word_length++] = stack[--stack_top];
      if (reducing)
        state[++valid_length] = s = s->transition(nm,g);
    }
    else
      break;
  }

  if (overflow)
  {
    if (failed)
    {
      String_Buffer sb;
      nm.rm.critical_status(2,reported ? 0 : 1,"Maximum word-length exceeded in reduction\n");
      *this = save;
      *failed = true;
      return false;
    }
    MAF_INTERNAL_ERROR(nm.maf.container,
                       ("Maximum word length exceeded. Was %d Now %d\n",
                       start_length,total_length));
    return false;
  }
  else
    set_length(word_length);
  if (failed && (time > 2 || start_length < nm.height() && nr_reductions > 50000))
    nm.rm.schedule_optimise(save);
  if (reported && reducing)
    nm.rm.critical_status(2,0,"Reduction completed.Time %lu. Reductions %lu."
                          " Length %d.\n",time,nr_reductions,word_length);
  return reducing;
}

/**/

Word_Length Equation_Word::slice(Word_Length start,Word_Length end)
{
  if (start || end < word_length)
  {
    if (end >= word_length)
      end = word_length;
    if (start >= end)
    {
      set_length(0);
      return 0;
    }

    if (start)
    {
      values += start;
      valid_length = 0;
      if (language == language_L1)
        language = language_L0;
      else if (language == language_L2)
        language = start >= 2 ? language_L0 : language_L1;
    }
    if (end != word_length)
    {
      if (language != language_A)
        language = language_L0;
    }
    word_length = end-start;
    if (valid_length > word_length)
      valid_length = word_length;
    values[word_length] = TERMINATOR_SYMBOL;
  }
#ifdef DEBUG
  format(&debug_word);
#endif
  return word_length;
}

/**/

Node_Reference Equation_Word::state_lhs(bool full_match) const
{
  Node_Reference answer;
  if (word_length < 2)
  {
    if (!word_length)
      return nm.start();
    answer = nm.start()->transition(nm,values[0]);
    state[1] = answer;
  }
  else
  {
    validate();
    if (!valid_length || valid_length+1 < word_length)
      return Node_Reference(0,0);
    answer = Node_Reference(nm,state[word_length]);
  }
  if (full_match && answer->length(nm) != word_length)
    answer = Node_Reference(0,0);
  return answer;
}

/**/

Node_Reference Equation_Word::state_word(bool full_match) const
{
  validate();
  if (valid_length < word_length)
    return Node_Reference(0,0);
  Node_Reference answer(nm,state[word_length]);
  if (full_match && answer->length(nm) != word_length)
    answer = Node_Reference(0,0);
  return answer;
}

/**/

Node_Reference Equation_Word::state_prefix() const
{
  validate();
  if (valid_length < word_length-1 || !word_length)
    return Node_Reference(0,0);
  Node_Reference answer(nm,state[word_length-1]);
  if (answer->length(nm)+1 != word_length)
    answer = Node_Reference(0,0);
  return answer;
}

/**/

Node_Reference Equation_Word::get_node(Word_Length length,bool lhs,
                                       Word_Length * new_prefix_length)
{
  Node_Reference nr = nm.start();

  Word_Length height = length + 1;
  Ordinal lvalue = values[0];
  if (nr->reduced.max_height < height)
    nr->node(nm).reduced.max_height = height;
  if (lhs)
    nr->node(nm).set_flags(NF_PREFIX_CHECKED);
  for (Word_Length i = 0; i < length;)
  {
    if (nr->is_final())
      MAF_INTERNAL_ERROR(nm.maf.container,("Attempting to create node with reducible prefix!\n"));
    Ordinal value = values[i++];
    Node_Reference child = nr->fast_child(nm,nr,value);
    if (child.is_null())
      child = nm.node_get()->node(nm).construct(nm,nr,lvalue,value,i);
    state[i] = nr = child;
    if (nr->reduced.max_height < height)
      nr->node(nm).reduced.max_height = height;

    if (new_prefix_length && !nr->flagged(NF_PREFIX_CHECKED))
    {
      *new_prefix_length = i;
      new_prefix_length = 0;
    }
    if (lhs)
      nr->node(nm).set_flags(NF_PREFIX_CHECKED);
  }
  if (nr->flagged(NF_NEW))
    nm.inspect(nr,false);

  valid_length = length;
  certificate = nm.current_certificate();
  return nr;
}

/**/

void Equation_Word::validate() const
{
  if (valid_length)
  {
    switch (nm.check_certificate(certificate))
    {
      case CV_VALID:
        break;
      case CV_INVALID:
        valid_length = 0;
        break;
      case CV_CHECKABLE:
        while (valid_length)
        {
          Word_Length v = Node_Reference(nm,state[valid_length])->valid_length(nm);
          if (v == valid_length)
            break;
          valid_length = Node_Reference(nm,state[valid_length])->length(nm) == valid_length ? v : valid_length-1;
        }
        break;
    }
  }
}

/**/

bool Equation_Word_Reducer::reduce(Word * word,const Word & start_word,
                                   unsigned flags,const FSA * )
{
  /* Convert word to reduced form by applying equations.
     Return value is 0 if word is reduced, or 1 if it gets reduced
  */
  Equation_Word::operator=(start_word);
  bool retcode = Equation_Word::reduce(flags & WR_PREFIX_ONLY ? AE_KEEP_LHS : 0);
  word->set_length(word_length);
  word->set_multiple(0,values,word_length);
  return retcode;
}

const Alphabet & Equation_Word::alphabet() const
{
  return nm.maf.alphabet;
}

