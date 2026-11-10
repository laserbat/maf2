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
$Log: equation.h $
Revision 1.9  2010/03/23 08:01:43Z  Alun
Changes for new style Node_Manager interface. BDT_Pool_Conjugate added,
and facility to discover if an equation was constructed by right conjugation
Revision 1.8  2009/09/12 18:48:25Z  Alun
printf() style functions cleaned up considerably
code now compiles cleanly in g++ both 32-bit and 64-bit
all tabs removed from source files


Revision 1.7  2008/12/24 18:41:54Z  Alun
Support for normal subgroup coset systems
Revision 1.6  2008/10/13 17:55:19Z  Alun
Early November checkin
Completely reworked pool and Difference_Tracker
Ported to Darwin
Revision 1.5  2008/08/25 06:32:56Z  Alun
"Early Sep 2008 snapshot"
Revision 1.3  2007/10/24 21:15:36Z  Alun
*/
#pragma once
#ifndef EQUATION_INCLUDED
#define EQUATION_INCLUDED 1
#ifndef MAFWORD_INCLUDED
#include "mafword.h"
#endif
#ifndef NODEBASE_INCLUDED
#include "nodebase.h"
#endif

//Classes referred too but defined elsewhere
class Container;
class Output_Stream;
class MAF;
class Presentation;
typedef Node Equation_Node;

struct Generator_Permutation
{
  Generator_Permutation * next;
  const Ordinal_Word * word;
  ~Generator_Permutation()
  {
    if (word)
      delete word;
  }
};

/**/

class Equation
{
  public:
    virtual ~Equation() {};
    virtual const Word & get_lhs() const = 0;
    virtual const Word & get_rhs() const = 0;
    void print(Container & container,Output_Stream * stream) const;
    // convert an equation to relator form. This will fail if the RHS is
    // not empty and not invertible
    bool relator(Ordinal_Word * relator,const Presentation & maf) const;
};

class Packed_Equation
{
  public:
    Packed_Word lhs_word;
    Packed_Word rhs_word;
    Packed_Equation(const Word &lhs_,const Word &rhs_) :
      lhs_word(lhs_),
      rhs_word(rhs_)
    {}
    Packed_Equation(const Equation & e) :
      lhs_word(e.get_lhs()),
      rhs_word(e.get_rhs())
   {}
};

class Linked_Packed_Equation : public Packed_Equation
{
  friend class Packed_Equation_List;
  private:
    Linked_Packed_Equation * next;
    ~Linked_Packed_Equation()
    {}
    Linked_Packed_Equation(const Word &lhs_,const Word &rhs_) :
      Packed_Equation(lhs_,rhs_),
      next(0)
    {}
    Linked_Packed_Equation(const Equation & e) :
      Packed_Equation(e),
      next(0)
    {}
  public:
    const Linked_Packed_Equation * get_next() const
    {
      return next;
    }
};

class Packed_Equation_List
{
  private:
    Linked_Packed_Equation *head;
    Linked_Packed_Equation *tail;
    unsigned nr_items;
  public:
    Packed_Equation_List() :
      head(0),
      tail(0),
      nr_items(0)
    {}
    ~Packed_Equation_List()
    {
      while (head)
      {
        Linked_Packed_Equation * next = head->next;
        delete head;
        head = next;
      }
    }
    void add(const Word &lhs,const Word &rhs)
    {
      Linked_Packed_Equation * new_lpe  = new Linked_Packed_Equation(lhs,rhs);
      if (!head)
        head = tail = new_lpe;
      else
      {
        tail->next = new_lpe;
        tail = new_lpe;
      }
      nr_items++;
    }
    const Linked_Packed_Equation * first() const
    {
      return head;
    }
    unsigned count() const
    {
      return nr_items;
    }
};

class Simple_Equation : public Equation
{
  public:
    Ordinal_Word lhs_word;
    Ordinal_Word rhs_word;
    Word_Length lhs_length;
    Word_Length rhs_length;
    Simple_Equation(const Node_Manager & nm,Equation_Handle e);
    Simple_Equation(const Word &lhs_,const Word &rhs_) :
      lhs_word(lhs_),
      rhs_word(rhs_)
    {
      lhs_length = lhs_word.length();
      rhs_length = rhs_word.length();
    }
    Simple_Equation(const Alphabet & alphabet,const Packed_Equation & pe) :
      lhs_word(alphabet,pe.lhs_word),
      rhs_word(alphabet,pe.rhs_word)
    {
      lhs_length = lhs_word.length();
      rhs_length = rhs_word.length();
    }
    Simple_Equation(const Alphabet &alphabet) :
      lhs_word(alphabet),
      rhs_word(alphabet),
      lhs_length(0),
      rhs_length(0)
    {}
    void read_equation(const Node_Manager &nm,Equation_Handle e);
    const Word & get_lhs() const
    {
      return lhs_word;
    }
    const Word & get_rhs() const
    {
      return rhs_word;
    }
    Total_Length total_length() const
    {
      return lhs_length + rhs_length;
    }
    Ordinal first_letter() const
    {
      return lhs_word.value(0);
    }
    Ordinal last_letter() const
    {
      return lhs_word.value(lhs_length-1);
    }
};

enum Basic_Derivation_Type
{
  BDT_Trivial,
  BDT_Axiom,
  BDT_Known,
  BDT_Pool,
  BDT_Right_Conjugate,
  BDT_Left_Conjugate,
  BDT_Conjunction,
  BDT_Equal_RHS,
  BDT_Reverse,
  BDT_Update,
  BDT_Half_Difference_1,
  BDT_Half_Difference_2,
  BDT_Inversion_Closure,
  BDT_Diff_Reduction,
  BDT_Double,
  BDT_Equal_LHS,
  BDT_Relator_Multiply,
  BDT_Normal_Closure,
  BDT_Pool_Conjugate,
  BDT_Unspecified
};

class Derivation
{
  protected:
    Basic_Derivation_Type type;
    union
    {
      struct
      {
        Node_ID e1;
        Node_ID e2;
        Word_Length offset;
      } overlap;
      struct
      {
        Node_ID e1;
        Ordinal g;
      } conjugate;
      struct
      {
        Node_ID state;
        Node_ID old_state;
        Transition_ID ti;
      } computation;
      ID pool_id;
    } data;
  public:
    /* Proof tracing only: the trace ids of equations a derivation cites,
       where they have to be captured when the derivation is created */
    int tpid[2];
    Derivation(Basic_Derivation_Type type_) :
      type(type_)
    {
      tpid[0] = tpid[1] = 0;
      data.overlap.e1 = 0;
      data.overlap.e2 = 0;
      data.overlap.offset = 0;
    }
    Derivation(Basic_Derivation_Type type_,Node_Handle e1) :
      type(type_)
    {
      tpid[0] = tpid[1] = 0;
      data.overlap.e1 = e1;
      data.overlap.e2 = 0;
      data.overlap.offset = 0;
    }
    Derivation(Basic_Derivation_Type type_,Node_Handle e1,
               Node_Handle e2,Word_Length offset = 0) :
      type(type_)
    {
      tpid[0] = tpid[1] = 0;
      data.overlap.e1 = e1;
      data.overlap.e2 = e2;
      data.overlap.offset = offset;
    }
    Derivation(Basic_Derivation_Type type_,Node_Handle e1,Ordinal g) :
      type(type_)
    {
      tpid[0] = tpid[1] = 0;
      data.conjugate.e1 = e1;
      data.conjugate.g = g;
    }
    Derivation(ID pool_id,Basic_Derivation_Type type_ = BDT_Pool) :
      type(type_)
    {
      tpid[0] = tpid[1] = 0;
      data.pool_id = pool_id;
    }
    Derivation(Basic_Derivation_Type type_,Node_Handle state,
               Transition_ID ti,
               Node_ID old_state = 0):
      type(type_)
    {
      tpid[0] = tpid[1] = 0;
      data.computation.state = state;
      data.computation.ti = ti;
      data.computation.old_state = old_state;
    }
    Node_ID right_conjugate() const
    {
      return type == BDT_Right_Conjugate ? data.conjugate.e1 : 0;
    }
    void set_trace(int pid1,int pid2)
    {
      tpid[0] = pid1;
      tpid[1] = pid2;
    }
    Basic_Derivation_Type get_type() const
    {
      return type;
    }

};
#endif
