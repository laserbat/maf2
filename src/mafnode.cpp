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


// $Log: mafnode.cpp $
// Revision 1.20  2010/06/10 13:57:33Z  Alun
// All tabs removed again
// Revision 1.19  2010/05/02 10:35:27Z  Alun
// Significant reworking required for new style Node_Manager which replaces
// Node * with Node_ID in Node members to save space in 64-bit MAF
// Revision 1.18  2009/11/09 21:29:59Z  Alun
// Minor changes to comments. Added dense_rm option
// Revision 1.17  2009/09/13 11:41:52Z  Alun
// printf() style functions cleaned up considerably
// code now compiles cleanly in g++ both 32-bit and 64-bit
// all tabs removed from source files
// Revision 1.16  2009/09/01 12:08:10Z  Alun_Williams
// Long comment updated to discuss sparse nodes
// Revision 1.15  2009/08/23 19:35:10Z  Alun
// stdio.h removed
// Revision 1.14  2009/01/01 21:20:40Z  Alun
// Now uses sparse transition table for nodes where this saves a significant amount of space
// Revision 1.14  2008/10/26 21:56:42Z  Alun
// Early November checkin
// Completely reworked pool and Difference_Tracker
// Ported to Darwin
// Revision 1.13  2008/10/01 06:57:08Z  Alun
// diagnostics removed
// Revision 1.12  2008/09/26 06:01:37Z  Alun
// Final version built using Indexer.
// Revision 1.7  2007/12/20 23:25:43Z  Alun
//

/*
   Although nodes are all share one C++ type, they represent several
   different types of data. In the original implementation of the code that
   became MAF there were three classes: State, Node, Equation, the latter two
   each being derived from the first, and considerable use was made of virtual
   functions. There is now only one class, partly because the four byte
   penalty of a vtable is unacceptable for a class which may require to be
   instantiated millions of times, but mainly because it is highly desirable
   that the information stored in a Node should not have to be moved when a
   word that was thought to be irreducible become reducible (the reverse
   case obviously never arises).

   Much of the node information is held in a union, so we need an
   indicator to tell us how to interpret it. The flags field is used for
   this. If neither NF_REDUCIBLE nor NF_REMOVED are present then the
   reduced information is present. If NF_REDUCIBLE is set then the
   reduction information is present, (or about to be set),
   but may be suspect if NF_REMOVED is also set. If NF_REMOVED is set on its
   own then the word is reducible, but no information on how to reduce it
   is present.

   If NF_REMOVED is set then the word has been found to have a reducible
   prefix and is gradually being destroyed.

   If a Node has forward links then these may either in one of two formats:

     "dense" format allows transitions to be computed quickly, because we
     store the next Node for every possible transition. In this case there
     are usually some links that are not links to a child, but to a trailing
     subword of the child word. These are known as prefix changes.

     "sparse" format is used when the number of generators is large, when
     storing the prefix change information would waste a lot of space,
     In this case only links to children are present and we have to look
     through the list of links to see if there is one for the child we
     are interested in. There is a single link to a "suffix" node which
     we must use to find other transitions. The "suffix" must contain a
     pointer to the longest trailing subword from which we may
     need to get a link to follow for a prefix change.

   Links are dense unless the NF_SPARSE flag is set. A Node is created
   NF_SPARSE if there are more than 8 possible children and it is not the
   root. However sparse links will later be turned into dense format
   if a Node acquires enough real children. When the NF_SPARSE flag is set
   determining the transition is slower because we first have to see if there
   is a link and if not look for the transition inside the suffix node, which
   may in turn need to look inside another suffix, and so on...

   Dense format is preferable because it makes the determination of
   transitions so much simpler. Since each node always has a pointer
   back to its prefix, it is easy to distinguish between a link to a
   "child" word and a "prefix change".
*/
#include "mafnode.h"
#include "maf_ew.h"
#include "equation.h"
#include "maf_nm.h"
#include "container.h"

/**/

Node_Reference Node::construct(Node_Manager & nm,
                               Node_Reference prefix,
                               Ordinal lvalue,
                               Ordinal rvalue_,
                               Word_Length word_length_)
{
#if MAF_USE_LOOKUP
  Node_ID id = reduction.node_id;
#else
  Node_ID id = this;
#endif
  Node_Reference answer(this,id);
  rvalue = rvalue_;
  reduced.word_length = word_length_;
  flags = NF_NEW;
  reduced.max_height = word_length_+1;
  lhs_node = 0;
  Ordinal child_end;
  Ordinal child_start;
  nm.valid_children(&child_start,&child_end,rvalue);
  if (nm.maf.options.dense_rm || child_end - child_start <= 8 ||
      word_length_ <= 1)
    reduced.child.dense = new Node_ID[child_end-child_start] - child_start;
  else
  {
    reduced.child.sparse = new Sparse_Node;
    flags |= NF_SPARSE;
  }

  reduced.inverse = 0;
  reduced.max_accepted = 0;
  reduced.lvalue = lvalue;
  if (word_length_)
  {
    this->prefix = prefix;
    /* When a node is first constructed it is not a prefix.
       Therefore the state to enter when the next character is
       seen should be the same as for the suffix. Links like this
       "forget" some initial part of the input, and have therefore
       accepted it. I call such links prefix changes. For word reduction
       to work correctly, there must be no prefix change links that
       cause a character to be accepted before all prefixes beginning
       at that character are ruled out. Thus if ABC is a node and BCD is
       a prefix of a LHS, but ABCD is not a prefix, the D link at ABC
       must point to BCD, and not to CD. If the equation involving BCD
       was not available when ABC was created, then the D link will initially
       be wrong, and must be corrected when the BCD LHS is inserted.
       That is one of the tasks performed in check_transitions(); */
    prefix->attach(nm)->set_flags(NF_PARENT);
    Node_Reference suffix = prefix->transition(nm,rvalue);
    MAF_ASSERT(!suffix->is_final(),nm.maf.container,("Suffix is reducible - error in prefix changes!\n"));
    if (prefix->flagged(NF_INVISIBLE) || suffix->flagged(NF_INVISIBLE))
      set_flags(NF_INVISIBLE);
    if (flagged(NF_SPARSE))
      reduced.child.sparse->suffix = suffix;
    else
    {
      /* Since this node was created dense its suffix must also be dense
         since either the suffix has the same possible children as this
         node or we are at _H */
      for (Ordinal g = child_start;g < child_end;g++)
      {
        Node_Reference child(nm,suffix->reduced.child.dense[g]);
        child->attach(nm)->set_flags(NF_TARGET);
        reduced.child.dense[g] = child;
        /* We no longer calculate reduced.max_accepted here because we always
           call inspect() on any new nodes created, and this calculates
           reduced.max_accepted anyway, and the calculation is now tricky */
        if (child->is_final())
          nm.publish_partial_reduction(Node_Reference(this,id),g);
      }
    }
    /* Attach the new node to its prefix */
    prefix->node(nm).link_replace(nm,prefix,rvalue,answer,suffix,LCF_Defer);
  }
  else
  {
    this->prefix = 0;
    attach(nm);
    for (Ordinal g = child_start;g < child_end;g++)
    {
      attach(nm);
      reduced.child.dense[g] = id;
    }
    reduced.max_accepted = 1;
    attach(nm);
    reduced.inverse = id;
  }
  nm.count_in(language_L0);
  return answer;
}

/**/

void Node::construct_equation(Node_Manager & nm,Node_ID nid,
                              Node_Reference prefix,Ordinal rvalue_,
                              unsigned flags_,Node_Reference rhs,
                              Node_Reference trailing_subword,bool new_node)
{
  /* This method is only called from Node_Manager::construct_equation()
     which already asserted that the prefix and rhs are not reducible */
  rhs->attach(nm);

  if (new_node)
  {
    this->prefix = prefix;
    prefix->attach(nm);
    rvalue = rvalue_;
  }
  else
  {
    if (reduction.reduced_node)
    {
      detach_rhs(nm,nid);
      if (reduction.trailing_subword)
        fast_find_node(nm,reduction.trailing_subword)->detach(nm,reduction.trailing_subword);
    }
    flags_ |= flags & (NF_TARGET|EQ_EXPANDED);
  }
  flags = flags_|NF_REDUCIBLE;
  reduction.queue_prev = 0;
  reduction.queue_next = 0;
  reduction.reduced_node = rhs;
  lhs_node = rhs->lhs_node;
  reduction.trailing_subword = trailing_subword;
  trailing_subword->attach(nm);
  rhs->node(nm).lhs_node = nid;
}

/**/

void Node::link_replace(Node_Manager & nm,Node_ID this_id,
                        Ordinal value,
                        Node_Reference new_child,
                        Node_Reference old_child,Link_Check_Flag check_flag)
{
  if (is_final())
    MAF_INTERNAL_ERROR(nm.maf.container,("Attempt to insert transition from reducible node!\n"));
  Ordinal child_start,child_end;
  nm.valid_children(&child_start,&child_end,rvalue);
  if (value < child_start || value >= child_end)
    MAF_INTERNAL_ERROR(nm.maf.container,("Attempt to insert invalid transition!\n"));
  if (new_child->length(nm) > reduced.word_length+1)
    MAF_INTERNAL_ERROR(nm.maf.container,("Attempt to insert state in wrong node!\n"));

  if (new_child != old_child)
  {
    if (value < 0 || new_child.is_null() ||
        value != new_child->rvalue && new_child != nm.start())
      MAF_INTERNAL_ERROR(nm.maf.container,("Attempt to insert invalid link %d %d\n",value,new_child->rvalue));
    if (new_child->prefix == this_id && old_child->prefix == this_id)
      MAF_INTERNAL_ERROR(nm.maf.container,("Attempt to insert duplicate node\n"));

    if (!flagged(NF_SPARSE))
    {
      new_child->attach(nm);
      if (old_child->prefix == this_id)
      {
        nm.revoke_transition(old_child->is_final());
        /* If we get here then we have detached the only real link into the
           child object, and it is now inaccessible, and we need to start the
           process of removing it altogether. But at this point there may be
           links into the object through prefix changes.
           If the call to link() comes from check_transitions() we can be
           confident that all these prefix changes links will get replaced by
           the time check_transitions completes. In any case the child will be
           in limbo until the last reference to it has gone, and if anything
           does try to the follow the child's transition links it will
           GPF/core dump */
        old_child->node(nm).remove_from_tree(nm,old_child);
      }
      old_child->detach(nm,old_child);
      reduced.child.dense[value] = new_child;
      if (new_child->prefix != this_id)
        new_child->node(nm).flags |= NF_TARGET;
    }
    else
    {
      Sparse_Node & sn = *reduced.child.sparse;
      if (new_child->prefix != this_id)
      {
        /* In this case we don't create a link, but we need to make sure
           our suffix is long enough to find the new child */
        if (new_child->parent(nm)->length(nm) > fast_find_node(nm,sn.suffix)->length(nm))
          sn.suffix = new_child->prefix;
      }
      else
      {
        new_child->attach(nm);

        Node_ID * new_children = sn.children;
        if (sn.nr_children == sn.nr_children_allocated)
          new_children = new Node_ID[++sn.nr_children_allocated];
        int i;
        for (i = sn.nr_children; i > 0;i--)
        {
          Node_ID child = sn.children[i-1];
          if (fast_find_node(nm,child)->rvalue < value)
            break;
          new_children[i] = child;
        }
        new_children[i] = new_child;
        if (new_children != sn.children)
        {
          while (--i >= 0)
            new_children[i] = sn.children[i];
          if (sn.children)
            delete [] sn.children;
          sn.children = new_children;
        }
        sn.nr_children++;
      }

      if (old_child->prefix == this_id)
      {
        int i;
        Node_ID old_id = old_child;
        for (i = 0; i < sn.nr_children ;i++)
          if (sn.children[i] == old_id)
            break;
        for (;i+1 < sn.nr_children;i++)
          sn.children[i] = sn.children[i+1];
        sn.nr_children--;
        nm.revoke_transition(old_child->is_final());
        old_child->node(nm).remove_from_tree(nm,old_child);
        old_child->detach(nm,old_child);
      }

      if (sn.nr_children * 3 > 2*(child_end - child_start))
      {
        /* This sparse node is rather full. Switch the node to being dense.
           We shall not bother with switching to sparse if children disappear
           later, as this is only likely to happen when a collapse is
           taking place */
        Node_ID * dense_child = new Node_ID[child_end-child_start] - child_start;
        for (Ordinal g = child_start;g < child_end;g++)
        {
          Node_Reference child = transition(nm,g);
          dense_child[g] = child;
          if (child->prefix != this_id)
            child->attach(nm)->set_flags(NF_TARGET);
        }
        delete reduced.child.sparse;
        reduced.child.dense = dense_child;
        clear_flags(NF_SPARSE);
      }
    }

    if (new_child->prefix == this_id)
      nm.publish_transition(new_child->is_final());
  }

  if (check_flag != LCF_Defer)
    nm.inspect(Node_Reference(this,this_id),check_flag == LCF_All);
}

/**/

void Node::rebind(Node * prefix_node,Node_Manager &nm,bool check_height) const
{
  if (!is_final())
  {
    prefix_node->set_flags(NF_PARENT);
    if (check_height && reduced.max_height > prefix_node->reduced.max_height)
      prefix_node->reduced.max_height = reduced.max_height;
    if (reduced.max_accepted > prefix_node->reduced.max_accepted)
      prefix_node->reduced.max_accepted = reduced.max_accepted;
  }
  else if (fast_is_primary())
  {
    /* although the input here is not really being accepted, we have to
       count it as though it were, since otherwise the trailing_subword
       might not be corrected in transition_check().
       The alternative would be to get rid of trailing_subword
       altogether and find it dynamically every time it is needed.
       The primary is said to have accepted
       (reduced.word_length+1)-trailing_subword->length()-1 letters:
       the +1 gets is to the length of the primaries LHS.
       the -1 is there because the trailing subword has to lose at least
       the first letter, so only letters after this can be said to be
       accepted.
    */
    Word_Length accepted = prefix_node->reduced.word_length -
                           fast_find_node(nm,reduction.trailing_subword)->length(nm);
    if (accepted > prefix_node->reduced.max_accepted)
      prefix_node->reduced.max_accepted = accepted;
  }
}

/**/

Node_Reference Node::fast_child(const Node_Manager &nm,Node_ID nid,Ordinal value) const
{
  if (flagged(NF_SPARSE))
    return sparse_child(nm,value);
  State answer = fast_find_node(nm,reduced.child.dense[value]);
  if (answer->prefix == nid)
    return Node_Reference(answer,reduced.child.dense[value]);
  return Node_Reference(0,0);
}

/**/

Node_Reference Node::sparse_child(const Node_Manager &nm,Ordinal value) const
{
  Sparse_Node & sn = *reduced.child.sparse;
  for (Ordinal i = 0; i < sn.nr_children;i++)
  {
    State child = fast_find_node(nm,sn.children[i]);
    if (child->rvalue >= value)
    {
      if (child->rvalue == value)
        return Node_Reference(child,sn.children[i]);
      break;
    }
  }
  return Node_Reference(0,0);
}

/**/

Node_Reference Node::sparse_transition(const Node_Manager &nm,Ordinal value) const
{
  Node_Reference answer = sparse_child(nm,value);
  if (answer.is_null())
  {
    State node = this;
    do
    {
      State suffix = fast_find_node(nm,node->reduced.child.sparse->suffix);
      if (suffix->flagged(NF_SPARSE))
      {
        answer = suffix->sparse_child(nm,value);
        node = suffix;
      }
      else
        return suffix->dense_transition(nm,value);
    }
    while (answer.is_null());
  }
  return answer;
}

/**/

void Node::make_reducible(Node_Manager & nm,Node_ID nid,unsigned short new_language)
{
  if (!flagged(NF_REDUCIBLE))
  {
    /* put all the equations that had this as the RHS into the list
       of equations whose RHS needs to be improved.
       We do this first because the RHS will appear to be valid at
       the moment*/
    for (Node_Reference lhs = first_lhs(nm);
         !lhs.is_null();
         lhs = lhs->next_lhs(nm))
      nm.publish_bad_rhs(lhs);

    /* We remove all our forward links. We have to be careful to do this
       otherwise inaccessible nodes would hang around in limbo for ever */
    if (!flagged(NF_SPARSE))
    {
      Ordinal child_start;
      Ordinal child_end;
      nm.valid_children(&child_start,&child_end,rvalue);
      for (Ordinal g = child_start;g < child_end;g++)
      {
        State child = fast_find_node(nm,reduced.child.dense[g]);
        if (child->prefix == nid)
        {
          nm.stats.ntc[child->is_final() ? 1 : 0]--;
          child->node(nm).remove_from_tree(nm,reduced.child.dense[g]);
        }
        child->detach(nm,reduced.child.dense[g]);
        reduced.child.dense[g] = 0;
      }
      delete [] (reduced.child.dense + child_start);
    }
    else
    {
      Sparse_Node & sn = *reduced.child.sparse;
      for (Ordinal i = 0; i < sn.nr_children;i++)
      {
        Node_Reference child(nm,sn.children[i]);
        nm.stats.ntc[child->is_final() ? 1 : 0]--;
        child->node(nm).remove_from_tree(nm,child);
        child->detach(nm,child);
        sn.children[i] = 0;
      }
      clear_flags(NF_SPARSE);
      delete &sn;
    }
    if (flagged(NF_IS_DIFFERENCE))
      nm.publish_bad_difference(Node_Reference(this,nid));
    Node_Reference inverse(nm,reduced.inverse);
    if (inverse)
    {
      if (!inverse->is_final())
        nm.publish_bad_inverse(inverse);
      inverse->detach(nm,inverse);
    }
    reduction.reduced_node = 0;
    reduction.id = 0;
    lhs_node = 0;
    reduction.trailing_subword = 0;
    reduction.expand_timestamp = 0;
    clear_flags(EQ_EXPANDED|EQ_PRIMARY);
    if (new_language != language_A)
    {
      nm.stats.ntc[0]--;
      nm.stats.ntc[1]++;
    }
  }
  nm.count_out(language(nm));
  if (new_language == language_A)
    nm.remove_node(Node_Reference(this,nid));
  else
  {
    flags |= NF_REDUCIBLE;
    if (new_language != language_L1)
      flags &= ~EQ_PRIMARY;
  }
  nm.count_in(new_language);
}

/**/

void Node::detach_rhs(Node_Manager & nm,Node_ID nid)
{
  /* Remove the equation from the list of equations with this RHS,
     if necessary. */
  Node_Reference rhs(nm,reduction.reduced_node);
  if (!rhs.is_null())
  {
    if (!rhs->is_final()) /* If the rhs is reducible we don't
                             care about the list */
    {
      if (rhs->lhs_node == nid)
        rhs->node(nm).lhs_node = lhs_node;
      else
      {
        const Equation_Node *e = find_node(nm,rhs->lhs_node);
        while (e && e->lhs_node != nid)
          e = find_node(nm,e->lhs_node);
        if (e)
          e->node(nm).lhs_node = lhs_node;
      }
    }
    rhs->detach(nm,reduction.reduced_node); /* Release our hold on the rhs */
  }
}

/**/

void Node::destroy(Node_Manager & nm,Node_ID nid)
{
  /* Destruction of nodes is complex, but is almost always caused by a prefix
     of the LHS becoming reducible.

     This function gets only gets called when there are no references left
     to a node, and we really are going to remove it.
  */
#ifdef DEBUG
  if (!flagged(NF_REMOVED))
    MAF_INTERNAL_ERROR(nm.maf.container,
                       ("Node::destroy() called for valid node %p\n",this));
#endif
  detach_rhs(nm,nid);
  if (prefix)
    fast_find_node(nm,prefix)->detach(nm,prefix); /* Release our hold on the prefix */
  nm.node_free(this);
}

/**/

int Node::compare(const Node_Manager & nm,State other) const
{
  if (this == other)
    return 0;
  const Alphabet & alphabet = nm.alphabet();
  Word_Length our_length = length(nm);
  Word_Length other_length = other->length(nm);
  if (alphabet.order_is_geodesic())
  {
    if (our_length > other_length)
      return 1;
    if (our_length < other_length)
      return -1;
    if (alphabet.order_is_shortlex())
    {
      Ordinal our_first = first_letter(nm);
      Ordinal other_first = other->first_letter(nm);
      if (our_first != other_first)
        return our_first > other_first ? 1 : -1;
    }
  }
  Ordinal_Word our_word(alphabet,our_length);
  Ordinal_Word other_word(alphabet,other_length);
  read_word(&our_word,nm);
  other->read_word(&other_word,nm);
  return our_word.compare(other_word);
}

/**/

int Node::compare(const Node_Manager &nm,const Word & other) const
{
  Word_Length our_length = length(nm);
  if (other.alphabet().order_is_geodesic())
  {
    Word_Length other_length = other.length();
    if (our_length > other_length)
      return 1;
    if (our_length < other_length)
      return -1;
    if (other.alphabet().order_is_shortlex())
    {
      Ordinal our_first = first_letter(nm);
      Ordinal other_first = other.value(0);
      if (our_first != other_first)
        return our_first > other_first ? 1 : -1;
    }
  }
  Ordinal_Word our_word(other.alphabet(),our_length);
  read_word(&our_word,nm);
  return our_word.compare(other);
}

/**/

void Node::print(Output_Stream * stream,const Node_Manager & nm,Node_ID nid) const
{
  Simple_Equation e(nm,Node_Reference(this,nid));
  e.print(nm.maf.container,stream);
}

/**/

Node_Reference Node::set_inverse(Node_Manager & nm,Node_ID nid,Node_Reference inverse,bool important)
{
  /* Although inverses should always be paired we do not rely on this,
     because it is perfectly possible that the inverse of the inverse
     of a word is not the original word while the Rewriter_Machine is
     not confluent. We set inverses as pairs when we are sure that we
     will get back to the same word, or when inversion is likely to be
     problematical, for example when recursive orderings are used. Otherwise
     we rely on the bad_inverse list to construct the inverse of the other
     word, and possibly a new equation as well.
     Inverses are not calculated unless word-differences are, or the
     appropriate tactics flag is set.
  */

  // When we have found the inverse of this word we will tell it to set its
  // own inverse if necessary
  if (!is_final())
  {
    Node_Reference self(this,nid);
    Node_Reference correct_inverse_inverse(self);
    Node_Reference old_inverse(nm,reduced.inverse);
    if (inverse.is_null())
      inverse = old_inverse;
    if (!inverse.is_null() && inverse->is_equation())
      inverse = inverse->fast_reduced_node(nm);
    if (inverse.is_null() || inverse->is_final())
    {
      Equation_Word iword(nm,(Word_Length) 0);
      if (read_inverse(&iword,nm))
      {
        bool failed = false;
        if (iword.reduce(0,&failed))
        {
          /* we can't be sure this is inverse's inverse */
          if (!nm.pd.inversion_difficult && !nm.is_confluent())
            correct_inverse_inverse = Node_Reference(0,0);
        }
        if (!failed)
        {
          inverse = iword.state_word();
          if (inverse.is_null() && (iword.length() < nm.height() || important))
            inverse = iword.get_node(iword.length(),false,0);
        }
        else
          inverse = Node_Reference(0,0);
      }
    }
    if (!inverse.is_null())
    {
      if (inverse != old_inverse)
      {
        inverse->attach(nm);
        if (!old_inverse.is_null())
        {
          nm.equate_nodes(old_inverse,inverse,
                          Derivation(BDT_Update,self,Transition_ID(-1)),-1,false,true);
          nm.publish_bad_inverse(old_inverse);
          old_inverse->detach(nm,old_inverse);
        }
        reduced.inverse = inverse;
      }
      Node_Reference inverse_inverse(nm,inverse->reduced.inverse);
      if (inverse_inverse != this)
      {
        if (inverse_inverse.is_null() || inverse_inverse->is_final())
        {
          inverse->node(nm).set_inverse(nm,inverse,correct_inverse_inverse,false);
          inverse_inverse = Node_Reference(nm,inverse->reduced.inverse);
          if (inverse_inverse.is_null())
          {
            /* If inverse still doesn't have an inverse, set it to
               this node. In this case reduction must have failed
               on the inverse of the inverse, or the creation of
               a node for a long word would have been required */
            inverse->node(nm).set_inverse(nm,inverse,self,false);
            inverse_inverse = Node_Reference(nm,inverse->reduced.inverse);
          }
        }

        if (inverse_inverse != this && !inverse_inverse.is_null() &&
            !nm.maf.options.no_equations)
        {
          nm.equate_nodes(self,inverse_inverse,
                          Derivation(BDT_Reverse,self,Transition_ID(-1)),
                          -1,false,true);
          /* our inverse already has an inverse, and it is not this node. So
             the two nodes are equal. It is not clear which belongs on the RHS
             (obviously we could find out, but it would involve reading and
             comparing the words for both nodes here, which is not convenient.
             If this is the node that belongs on the LHS it will detach itself
             from the inverse automatically when it is made reducible.
             But, if it belongs on the RHS our inverse won't know how to
             tell us to detach from it if it becomes reducible. So we have to
             set the flag that says we are a bad inverse. update_machine() will
             create the equation above before checking our inverse, so it will
             be able to correct things. */
          nm.publish_bad_inverse(self);
        }
      }
    }
    else if (!old_inverse.is_null())
    {
      nm.publish_bad_inverse(old_inverse);
      old_inverse->detach(nm,old_inverse);
      reduced.inverse = 0;
    }
    return inverse;
  }
  return Node_Reference(0,0);
}

/**/

Node_Reference Node::inverse(Node_Manager &nm,Node_ID nid) const
{
  Node_Reference answer = reduced_node(nm,nid);
  if (!answer->is_final() && answer->reduced.inverse)
  {
    answer = Node_Reference(nm,answer->reduced.inverse);
    if (!answer->is_final())
      return answer;
    answer = answer->reduced_node(nm,answer);
    if (!answer->is_final())
      return answer;
  }
  return Node_Reference(0,0);
}

/**/

Node_Reference Node::left_half_difference(Node_Manager & nm,Node_ID nid,Ordinal value) const
{
  /* invert and multiply on right. Answer is therefore inverse
     of Xw. Called twice with x,y we get Xwy */
  Node_Reference answer = inverse(nm,nid);
  if (!answer.is_null() && value != PADDING_SYMBOL)
  {
    if (answer->prefix && value == nm.inverse(answer->rvalue))
      answer = answer->parent(nm);
    else
    {
      answer = answer->fast_child(nm,answer,value);
      if (!answer.is_null() && answer->is_equation())
        answer = answer->fast_reduced_node(nm);
    }
  }
  return answer;
}

/**/

Node_Reference Node::right_half_difference(Node_Manager &nm,Node_ID nid,Ordinal value) const
{
  /* right multiply and invert. Answer is therefore inverse of wy.
     Called twice we get Xwy.*/
  Node_Reference answer = value==-1 ? Node_Reference(nm,nid) :
                             (prefix && value==nm.inverse(rvalue) ?
                              parent(nm) : child(nm,nid,value));
  if (!answer.is_null())
    answer = answer->inverse(nm,answer);
  return answer;
}

/**/

Node_Reference Node::difference(Node_Manager & nm,Node_ID nid,Ordinal lvalue,Ordinal rvalue) const
{
  Node_Reference answer = left_half_difference(nm,nid,lvalue);
  if (!answer.is_null())
  {
    answer = answer->left_half_difference(nm,answer,rvalue);
    if (!answer.is_null())
      return answer;
  }

  answer = right_half_difference(nm,nid,rvalue);
  if (!answer.is_null())
    answer = answer->right_half_difference(nm,answer,lvalue);
  return answer;
}

bool Node::read_inverse_values(Ordinal * values,const Node_Manager & nm) const
{
  Word_Length i = 0;
  State node = this;
  Word_Length length = node->length(nm);
  while (length--)
  {
    if ((values[i++] = nm.inverse(node->rvalue))==INVALID_SYMBOL)
      return false;
    node = find_node(nm,node->prefix);
  }
  return true;
}

/**/

Language Node::language(const Node_Manager & nm) const
{
  if (!is_final())
    return language_L0;
  if (flagged(NF_REMOVED))
    return language_A;
  if (flagged(EQ_PRIMARY))
    return language_L1;
  /* The prefix is currently irreducible since otherwise we would
     have been removed already.
     The trailing subword is reducible, but we cannot
     guarantee it has an irreducible prefix here, so we have to
     find its length with length() rather than lhs_length() */
  if (fast_find_node(nm,reduction.trailing_subword)->length(nm) ==
      fast_find_node(nm,prefix)->reduced.word_length)
    return language_L2;
  return language_L3;
}

/**/

Node_Status Node::status() const
{
  if (!flagged(NF_REDUCIBLE))
    return flagged(NF_REMOVED) ? NS_Removed_Node : NS_Irreducible;
  if (flagged(NF_REMOVED))
    return NS_Removed_Equation;
  switch (flags & EQ_State_Mask)
  {
    case EQ_Unconjugated:
      return NS_Unconjugated_Equation;
    case EQ_Undifferenced:
      return NS_Undifferenced_Equation;
    case EQ_Dubious_Uncorrected:
      return NS_Dubious_Uncorrected_Equation;
    case EQ_Dubious_Unconjugated:
      return NS_Dubious_Unconjugated_Equation;
    case EQ_Dubious_Undifferenced:
      return NS_Dubious_Undifferenced_Equation;
    case EQ_Dubious_Secondary:
      return NS_Dubious_Secondary_Equation;
    case EQ_Uncorrected:
      return NS_Uncorrected_Equation;
    case EQ_Repair_Pending:
      return NS_Correction_Pending_Equation;
    case EQ_Weak_Secondary:
      return NS_Weak_Secondary_Equation;
    case EQ_Adopted:
      return NS_Adopted_Equation;
    case EQ_Expanded:
      return NS_Expanded_Equation;
    case EQ_Oversized:
      return NS_Oversized_Equation;
    default:
      return NS_Invalid;
  }
}

/**/

bool Node::set_status(Node_Manager &nm,Node_Status new_status)
{
  Node_Status was_status = status();
  if (!flagged(NF_REDUCIBLE))
  {
    if (new_status == NS_Removed_Node)
    {
      flags |= NF_REMOVED;
      return true;
    }
    else
    {
      MAF_INTERNAL_ERROR(nm.maf.container,
                         ("Equation being given unexpected status %d in"
                          " Node::set_status()\n",new_status));
      return false;
    }
  }
  ID new_id = id();
  flags &= ~EQ_State_Mask;
  switch (new_status)
  {
    case NS_Expanded_Equation:
      if (flags & EQ_QUEUED)
        MAF_INTERNAL_ERROR(nm.maf.container,
                           ("Cannot make queued equation with status %d"
                            " expandable in Node::set_status()\n",was_status));

      flags |= EQ_Expanded;
      if (nm.maf.options.strategy & MSF_SELECTIVE_PROBE)
        new_id = 0;
    case NS_Adopted_Equation:
      reduction.expand_timestamp = 0;
      reduction.id = new_id ? new_id : ++nm.last_id;
      flags |= EQ_Adopted;
      return true;
    case NS_Weak_Secondary_Equation:
      flags |= EQ_Weak_Secondary;
      return true;
    case NS_Undifferenced_Equation:
      flags |= EQ_Undifferenced;
      return true;
    case NS_Unconjugated_Equation:
      flags |= EQ_Unconjugated;
      return true;
    case NS_Dubious_Unconjugated_Equation:
      flags |= EQ_Dubious_Unconjugated;
      return true;
    case NS_Dubious_Undifferenced_Equation:
      flags |= EQ_Dubious_Undifferenced;
      return true;
    case NS_Dubious_Secondary_Equation:
      flags |= EQ_Dubious_Secondary;
      return true;
    case NS_Uncorrected_Equation:
      flags |= EQ_Uncorrected;
      return true;
    case NS_Dubious_Uncorrected_Equation:
      flags |= EQ_Dubious_Uncorrected;
      return true;
    case NS_Correction_Pending_Equation:
      flags |= EQ_Repair_Pending;
      return true;
    case NS_Removed_Equation:
      flags = (flags & ~EQ_PRIMARY) | EQ_Removed | NF_REMOVED;
      reduction.removed_id = ++nm.last_removed_id;
      return true;
    case NS_Oversized_Equation:
      flags |= EQ_Oversized;
      return true;

    case NS_Invalid:      // cases handled explicitly to tell fussy compilers we
    case NS_Removed_Node: // have thought about them
    case NS_Irreducible:
    default:
      MAF_INTERNAL_ERROR(nm.maf.container,
                         ("Equation being given unexpected status"
                          " %d in Node::set_status\n",new_status));
      return false;
  }
}

/**/

Container & Node::container(const Node_Manager &nm)
{
  return nm.maf.container;
}

/**/

Node_Iterator::Node_Iterator(const Node_Manager & nm_,
                             //Node_Filter * filter_,
                             Node_Reference root_) :
  nm(nm_),
  current_word(nm_.alphabet(),nm_.height())
{
  begin(root_);
}

/**/

void Node_Iterator::begin(Node_Reference root_)
{
  root = root_.is_null() ? nm.start() : root_;
  current = root;
  position = 1;
  parent = root->parent(nm);
  depth = root->length(nm);
  root->read_word(&current_word,nm);
}

/**/

bool Node_Iterator::next(Node_Reference *nr,bool inside)
{
  if (!inside)
  {
    if (current == root)
    {
      *nr = current = Node_Reference(0,0);
      position = 3;
      return false;
    }
    current = parent;
    parent = parent->parent(nm);
    depth--;
    position = 0;
  }

  for (;;)
  {
    parent = current;
    if (!parent->is_final())
    {
      if (!(parent->flags & NF_SPARSE))
      {
        nm.valid_children(&child_start,&child_end,parent->last_letter());
        Ordinal g = position!=0 ? child_start : current_word.value(depth)+1;
        for (;g < child_end;g++)
        {
          current = parent->dense_transition(nm,g);
          if (current->prefix == parent)
          {
            current_word.set_length(depth+1);
            current_word.set_code(depth++,g);
            *nr = current;
            position = 2;
            return true;
          }
        }
      }
      else
      {
        Node::Sparse_Node &sn = *parent->reduced.child.sparse;
        if (position != 0)
        {
          if (sn.nr_children)
          {
            current = Node_Reference(nm,sn.children[0]);
            current_word.set_length(depth+1);
            current_word.set_code(depth++,current->last_letter());
            *nr = current;
            position = 2;
            return true;
          }
        }
        else
        {
          Ordinal g = current_word.value(depth);
          for (Ordinal i = 0; i < sn.nr_children;i++)
          {
            current = Node_Reference(nm,sn.children[i]);
            if (current->last_letter() > g)
            {
              current_word.set_length(depth+1);
              current_word.set_code(depth++,current->last_letter());
              *nr = current;
              position = 2;
              return true;
            }
          }
        }
      }
    }
    if (parent == root)
    {
      *nr = current = Node_Reference(0,0);
      position = 3;
      return false;
    }
    current = parent->parent(nm);
    depth--;
    position = 0;
  }
}

/**/

struct Candidate
{
  Node_Reference e2;
  Total_Length l;
  Total_Length r;
  Total_Length tl;
  Total_Length overlap_length;
  bool seen_before;
  ID id;
};

void Node::find_candidates(const Find_Candidate & find,Node_ID nid,Total_Length overlap_length) const
{
  /* The process of creating an equation can result in major changes
     to the structure of the Node_Manager tree. However, because
     nodes for equations that have been removed from the tree are
     not actually destroyed during a single expansion cycle we can be
     sure that any equations do still exist while we are processing a
     list of candidates.
     Therefore we can safely build a list of plausible pairs and then
     process it one by one in expand().
  */
  Word_Length word_length = this->length(find.nm);
  Word_Length offset = overlap_length - word_length;
  /* The length of the node is only relevant for checking that
     prefix changes have not taken us past the end of the original LHS,
     and must not be used to exit from the recursion because its length
     exceeds the length of equations we are trying to pair.
     It is perfectly possible for node length to exceed the length of the
     equations we are looking to pair at the moment, because a prefix change
     link from the node could take us to a shorter equation. */

  if (offset < find.lhs_length)
  {
    if (is_final())
    {
      if (is_adopted())
      {
        Candidate c;

        c.l = find.rhs_length + overlap_length - find.lhs_length;
        c.r = reduced_length(find.nm) + offset;
        c.tl = c.l+c.r;
        c.id = id();
        c.seen_before = overlap_length <= find.upper_overlap_limit &&
                        c.id < find.min_id &&
                        find.strategy & MSF_RESPECT_RIGHT &&
                        (overlap_length <= find.old_lower_overlap_limit ||
                         c.tl <= find.old_check_limit);

        c.overlap_length = overlap_length;
        c.e2 = Node_Reference(this,nid);
        bool already_present = c.seen_before &&
                               c.tl <= find.old_keep_limit &&
                               overlap_length <= find.old_lower_overlap_limit;
        if (!already_present)
        {
          bool will_see_again = c.tl > find.keep_limit ||
                                find.strategy & MSF_SELECTIVE_PROBE &&
                                status() != NS_Expanded_Equation;
          if (c.id <= find.max_id && !will_see_again ||
              find.is_interesting(c))
          {
            attach(find.nm);
            find.candidates_list->add_state(c.e2,offset);
          }
          else
            find.nm.overlap_filter.set_filtered();
        }
      }
      if (!fast_is_primary())
      {
        Node_ID child = reduction.trailing_subword;
        fast_find_node(find.nm,child)->find_candidates(find,child,
                                                       overlap_length);
      }
      else if (find.strategy & MSF_DEEP_RIGHT)
      {
        /* generally we only look at overlaps up to the first overlapping
           reduction. But optionally we can consider all overlaps */
        fast_find_node(find.nm,reduction.trailing_subword)->find_candidates(find,reduction.trailing_subword,overlap_length);
      }
    }
    else
    {
      if (overlap_length < find.upper_overlap_limit)
      {
        Ordinal child_end;
        Ordinal child_start;
        find.nm.valid_children(&child_start,&child_end,rvalue);
        if (overlap_length - word_length + 1 < find.lhs_length)
        {
          if (!flagged(NF_SPARSE))
            for (Ordinal g = child_start; g < child_end;g++)
            {
              Node_Reference child = dense_transition(find.nm,g);
              child->find_candidates(find,child,overlap_length+1);
            }
          else
            for (Ordinal g = child_start; g < child_end;g++)
            {
              Node_Reference child = sparse_transition(find.nm,g);
              child->find_candidates(find,child,overlap_length+1);
            }
        }
        else
        {
          for (Ordinal g = child_start; g < child_end;g++)
          {
            Node_Reference child = fast_child(find.nm,nid,g);
            if (child)
              child->find_candidates(find,child,overlap_length+1);
          }
        }
      }
//      else
//      find.nm.overlap_filter.set_filtered();
    }
  }
}

/**/

Find_Candidate::Find_Candidate(Node_Manager & nm_,
                               Node_List * candidates_list_,
                               ID min_id_,
                               ID max_id_,
                               Equation_Handle e1_, Word_Length suffix_length,
                               const Equation_Word & e1_word_,
                               Overlap_Search_Type ost_,
                               bool favour_differences_) :
  nm(nm_),
  candidates_list(candidates_list_),
  min_id(min_id_),
  max_id(max_id_),
  lhs_length(e1_word_.length()),
  rhs_length(e1_->reduced_length(nm)),
  e1(e1_),
  word(e1_word_),
  ost(ost_),
  favour_differences(favour_differences_)
{
  if (ost <= OST_Minimal)
    suffix_length = nm.overlap_filter.probe_depth;
  lower_overlap_limit = upper_overlap_limit = lhs_length + suffix_length;
  probe_style = nm.maf.options.probe_style;
  if (!(nm.maf.options.filters & 1))
    lower_overlap_limit = old_lower_overlap_limit = upper_overlap_limit = old_upper_overlap_limit = UNLIMITED;
  else
  {
    if (lower_overlap_limit > nm.overlap_filter.overlap_limit ||
        probe_style & 2)
      lower_overlap_limit = nm.overlap_filter.overlap_limit;
    if (probe_style & 2 && ost > OST_Minimal ||
        upper_overlap_limit >= nm.overlap_filter.forbidden_limit)
      upper_overlap_limit = nm.overlap_filter.forbidden_limit-1;
    old_lower_overlap_limit = old_upper_overlap_limit = lhs_length + nm.old_overlap_filter.probe_depth;
    if (old_lower_overlap_limit > nm.old_overlap_filter.overlap_limit ||
        probe_style & 2)
      old_lower_overlap_limit = nm.old_overlap_filter.overlap_limit;
    /* The previous pass searched no further than one short of its
       forbidden limit (see upper_overlap_limit above), even where its
       overlap_limit was as long: an overlap it never looked at must not
       be taken as considered. */
    if (nm.old_overlap_filter.forbidden_limit &&
        old_lower_overlap_limit >= nm.old_overlap_filter.forbidden_limit)
      old_lower_overlap_limit = nm.old_overlap_filter.forbidden_limit-1;
    if (probe_style & 2 ||
        old_upper_overlap_limit >= nm.old_overlap_filter.forbidden_limit)
      old_upper_overlap_limit = nm.old_overlap_filter.forbidden_limit-1;
  }
  keep_limit = nm.overlap_filter.keep_limit;
  old_keep_limit = nm.old_overlap_filter.keep_limit;
  vital_limit = nm.overlap_filter.vital_limit;
  if (nm.maf.options.filters & 2)
  {
    check_limit = nm.overlap_filter.check_limit;
    old_check_limit = nm.old_overlap_filter.check_limit;
    if (!(nm.maf.options.filters & 4))
      lower_overlap_limit = 0;
  }
  else
  {
    check_limit = old_check_limit = UNLIMITED;
    if (probe_style & 1)
    {
      lower_overlap_limit = upper_overlap_limit;
      old_lower_overlap_limit = old_upper_overlap_limit;
    }
    else
    {
      upper_overlap_limit = lower_overlap_limit;
      old_upper_overlap_limit = old_lower_overlap_limit;
    }
  }

  if (favour_differences)
  {
    old_check_limit = old_keep_limit;
    if (ost != OST_Full)
      check_limit = 0;
  }

  strategy = nm.maf.options.strategy;
  e1_special = 0;
  if (favour_differences)
  {
    if (e1->flagged(EQ_INTERESTING))
      e1_special = 1;
    if (e1->flagged(EQ_HAS_PRIMARY_DIFFERENCES))
      e1_special |= 2;
  }
  if (ost != OST_Full)
    strategy |= MSF_RESPECT_RIGHT;
}

/**/

bool Find_Candidate::is_interesting(Candidate &c) const
{
  if (strategy & MSF_USE_ERAS && c.id > max_id)
    return false;
  if (c.tl <= check_limit)
    return true;
  if (c.overlap_length <= lower_overlap_limit)
    return true;
  if (!favour_differences || c.seen_before)
    return false;
  switch (ost)
  {
    case OST_Aborted_Pass:
    case OST_Minimal:
      return false;
    case OST_Quick:
      return e1_special &&
             c.e2->flagged(EQ_HAS_PRIMARY_DIFFERENCES+EQ_INTERESTING) &&
             (e1_special >= 2 || c.e2->flagged(EQ_HAS_PRIMARY_DIFFERENCES));
    case OST_Moderate:
      return e1_special && c.e2->flagged(EQ_INTERESTING+EQ_HAS_PRIMARY_DIFFERENCES);
    case OST_Full:
      return e1_special >= 2 ||
             c.e2->flagged(EQ_HAS_PRIMARY_DIFFERENCES) ||
             e1_special && c.e2->flagged(EQ_INTERESTING);
  }
  return false;
}
