/*
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

/* Proof tracing - see maf_trace.h for what the trace contains.

   The work here is turning what MAF records about an equation into
   rewriting steps.  An equation starts as two words known to be equal,
   with a justification that depends on where it came from - the two
   reductions of an overlap, say - and is then rewritten, swapped, and
   in groups and other cancellative monoids also cancelled and balanced.
   Full_Derivation records each of those operations as it happens (see
   maf_we.h), and emit_equation() replays them on the words, keeping a
   chain of steps from the current left hand side to the current right
   hand side:

     reduce the lhs by p        prepend p, backwards
     reduce the rhs by p        append p
     swap the sides             reverse the chain
     cancel a common suffix s   insert s*w at the end of the lhs, where
                                s*w = 1 is an equation we know, run the old
                                chain, and cancel s*w again; likewise with
                                w*s = 1 for a common prefix
     balance, unbalance         the same, with g*inv(g) = 1

   which is the same arithmetic KBMAG's kbtrace.c does for its insert().
*/

#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include "awcc.h"
#include "maf.h"
#include "maf_rm.h"
#include "mafnode.h"
#include "maf_nm.h"
#include "maf_we.h"
#include "hash.h"
#include "mafword.h"
#include "maf_rws.h"
#include "container.h"
#include "maf_trace.h"

Proof_Trace * proof_trace = 0;
static void read_word(Trace_Word * tw,const Word & w);
static std::string trace_filename;

/**/

void Proof_Trace::set_filename(const char * filename)
{
  trace_filename = filename ? filename : "";
}

/**/

const char * Proof_Trace::wanted()
{
  return trace_filename.empty() ? 0 : trace_filename.c_str();
}

/**/

Proof_Trace::Proof_Trace(const char * filename_) :
  file(0),
  filename(filename_),
  next_pid(1),
  gaps(0),
  debug(getenv("MAF_TRACE_DEBUG") != 0),
  nr_generators(0),
  have_final(false),
  final_confluent(false),
  final_lostinfo(false),
  resumed(false),
  closed(false)
{
  lhs_of.push_back(Trace_Word());
  rhs_of.push_back(Trace_Word());
}

/**/

Proof_Trace::~Proof_Trace()
{
  close();
}

/**/

static int kbmag_ordering(Word_Ordering wo)
{
  /* KBMAG's numbering - SHORTLEX, RECURSIVE, RT_RECURSIVE, WTLEX,
     WREATHPROD, NONE - which is what the trace format records.  MAF
     orderings KBMAG does not have come out as NONE */
  switch (wo)
  {
    case WO_Shortlex:
      return 0;
    case WO_Recursive:
      return 1;
    case WO_Right_Recursive:
      return 2;
    case WO_Weighted_Shortlex:
      return 3;
    case WO_Wreath_Product:
      return 4;
    default:
      return 5;
  }
}

/**/

bool Proof_Trace::start(MAF & maf)
{
  if (maf.properties().is_coset_system)
  {
    fprintf(stderr,"Proof tracing is not supported for coset systems.\n");
    return false;
  }
  if (maf.has_permutations())
  {
    /* insert_axiom() would add the image of each axiom under the
       permutations as further axioms, which the presentation does not
       state */
    fprintf(stderr,"Proof tracing is not supported with generator"
                   " permutations.\n");
    return false;
  }
  file = fopen(filename.c_str(),"w");
  if (!file)
  {
    fprintf(stderr,"Error: cannot open trace file %s\n",filename.c_str());
    return false;
  }
  const Alphabet & alphabet = maf.alphabet;
  nr_generators = alphabet.letter_count();
  inverse_of.resize(nr_generators);
  right_inverse.assign(nr_generators,0);
  left_inverse.assign(nr_generators,0);
  fprintf(file,"#KBTRACE 1\n");
  fprintf(file,"#ordering %d\n",kbmag_ordering(alphabet.order_type()));
  fprintf(file,"#resumed %d\n",resumed ? 1 : 0);
  fprintf(file,"#gens %d",int(nr_generators));
  for (Ordinal g = 0; g < nr_generators;g++)
    fprintf(file," %s",alphabet.glyph(g).string());
  fputc('\n',file);
  fprintf(file,"#inverses %d",int(nr_generators));
  for (Ordinal g = 0; g < nr_generators;g++)
  {
    inverse_of[g] = maf.inverse(g);
    fprintf(file," %d",inverse_of[g] == INVALID_SYMBOL ? 0 : inverse_of[g]+1);
  }
  fputc('\n',file);
  /* What the trace may take as given: the axioms of the presentation,
     which given_axiom() is told about as they are read, and g*inv(g) = 1 */
  for (Ordinal g = 0; g < nr_generators;g++)
    if (inverse_of[g] != INVALID_SYMBOL)
    {
      Trace_Word gG(2);
      gG[0] = g;
      gG[1] = inverse_of[g];
      given.insert(key(gG,Trace_Word()));
    }
  return true;
}

/**/

void Proof_Trace::given_axiom(const Word & lhs,const Word & rhs)
{
  Trace_Word l,r;
  read_word(&l,lhs);
  read_word(&r,rhs);
  given.insert(key(l,r));
  given.insert(key(r,l));
}

/**/

void Proof_Trace::close()
{
  if (!file || closed)
    return;
  closed = true;
  fprintf(file,"#confluent %d\n",final_confluent ? 1 : 0);
  fprintf(file,"#lostinfo %d\n",final_lostinfo ? 1 : 0);
  if (have_final)
  {
    fprintf(file,"F %d",int(final_pids.size()));
    for (size_t i = 0; i < final_pids.size();i++)
      fprintf(file," %d",final_pids[i]);
    fputc('\n',file);
  }
  fclose(file);
  file = 0;
  if (gaps)
    fprintf(stderr,"#Warning: %lu step(s) of the trace could not be justified;"
                   " the proof will be incomplete.\n",gaps);
}

/**/

static void read_word(Trace_Word * tw,const Word & w)
{
  Word_Length l = w.length();
  const Ordinal * values = w.buffer();
  tw->assign(values,values+l);
}

/**/

static void read_node(Trace_Word * tw,const Node_Manager & nm,Node_ID n)
{
  Ordinal_Word ow(nm.maf.alphabet);
  Node::fast_find_node(nm,n)->read_word(&ow,nm);
  read_word(tw,ow);
}

/**/

std::string Proof_Trace::key(const Trace_Word & lhs,const Trace_Word & rhs)
{
  std::string answer;
  answer.reserve((lhs.size()+rhs.size()+1)*sizeof(Ordinal));
  if (!lhs.empty())
    answer.append((const char *) &lhs[0],lhs.size()*sizeof(Ordinal));
  Ordinal separator = -1;
  answer.append((const char *) &separator,sizeof(Ordinal));
  if (!rhs.empty())
    answer.append((const char *) &rhs[0],rhs.size()*sizeof(Ordinal));
  return answer;
}

/**/

int Proof_Trace::content(const Trace_Word & lhs,const Trace_Word & rhs) const
{
  std::unordered_map<std::string,int>::const_iterator i = by_content.find(key(lhs,rhs));
  return i == by_content.end() ? 0 : i->second;
}

/**/

int Proof_Trace::equation(const Word & lhs,const Word & rhs) const
{
  Trace_Word l,r;
  read_word(&l,lhs);
  read_word(&r,rhs);
  return content(l,r);
}

/**/

int Proof_Trace::node_equation(const Node_Manager & nm,Node_ID n)
{
  State node = Node::fast_find_node(nm,n);
  if (!node->is_equation())
    return 0;
  Trace_Word lhs,rhs;
  read_node(&lhs,nm,n);
  read_node(&rhs,nm,node->raw_reduced_node(nm));
  std::unordered_map<Node_ID,int>::const_iterator i = node_pid.find(n);
  if (i != node_pid.end() && lhs_of[i->second] == lhs && rhs_of[i->second] == rhs)
    return i->second;
  return content(lhs,rhs);
}

/**/

void Proof_Trace::node_chain(const Node_Manager & nm,Node_ID n,int pos,
                             Trace_Chain * chain)
{
  /* A reduction by n replaces lhs(n) by fast_reduced_node(), which follows
     the raw right hand sides through any that are themselves left hand
     sides; each of those is a rewrite of its own at the same place */
  for (;;)
  {
    Trace_Step step;
    step.pid = node_equation(nm,n);
    step.pos = pos;
    step.dir = 1;
    chain->push_back(step);
    State node = Node::fast_find_node(nm,n);
    if (!node->is_equation())
      break;
    Node_Reference next = node->raw_reduced_node(nm);
    if (!next->is_equation())
      break;
    n = next;
  }
}

/**/

void Proof_Trace::put_word(const Ordinal * values,size_t length)
{
  if (!length)
  {
    fputc('.',file);
    return;
  }
  fprintf(file,"%d",values[0]+1);
  for (size_t i = 1; i < length;i++)
    fprintf(file,",%d",values[i]+1);
}

/**/

int Proof_Trace::write(const char * origin,const Trace_Word & lhs,
                       const Trace_Word & rhs,const Trace_Chain & chain)
{
  int pid = next_pid++;
  fprintf(file,"E %d ",pid);
  put_word(lhs);
  fputc(' ',file);
  put_word(rhs);
  fprintf(file," %s :",origin);
  for (size_t i = 0; i < chain.size();i++)
  {
    fprintf(file," %d@%d%c",chain[i].pid,chain[i].pos,chain[i].dir > 0 ? '>' : '<');
    if (!chain[i].pid)
      gaps++;
  }
  fputc('\n',file);
  lhs_of.push_back(lhs);
  rhs_of.push_back(rhs);
  by_content[key(lhs,rhs)] = pid;
  if (rhs.empty() && !lhs.empty())
  {
    right_inverse[lhs.front()] = pid;
    left_inverse[lhs.back()] = pid;
  }
  return pid;
}

/**/

bool Proof_Trace::apply(Trace_Word * w,const Trace_Step & step,std::string * why) const
{
  if (step.pid <= 0 || step.pid >= next_pid)
  {
    *why = "no equation";
    return false;
  }
  const Trace_Word & from = step.dir > 0 ? lhs_of[step.pid] : rhs_of[step.pid];
  const Trace_Word & to = step.dir > 0 ? rhs_of[step.pid] : lhs_of[step.pid];
  if (step.pos < 0 || step.pos + from.size() > w->size() ||
      !std::equal(from.begin(),from.end(),w->begin()+step.pos))
  {
    char buffer[64];
    sprintf(buffer,"step %d@%d%c does not match",step.pid,step.pos,
            step.dir > 0 ? '>' : '<');
    *why = buffer;
    return false;
  }
  w->erase(w->begin()+step.pos,w->begin()+step.pos+from.size());
  w->insert(w->begin()+step.pos,to.begin(),to.end());
  return true;
}

/**/

bool Proof_Trace::replay(Trace_Word w,const Trace_Chain & chain,
                         const Trace_Word & target,std::string * why) const
{
  for (size_t i = 0; i < chain.size();i++)
    if (!apply(&w,chain[i],why))
      return false;
  if (w != target)
  {
    *why = "derivation does not reach the right hand side";
    return false;
  }
  return true;
}

/**/

static void reverse_chain(Trace_Chain * chain)
{
  Trace_Chain answer(chain->rbegin(),chain->rend());
  for (size_t i = 0; i < answer.size();i++)
    answer[i].dir = -answer[i].dir;
  chain->swap(answer);
}

/**/

static void push_step(Trace_Chain * chain,int pid,int pos,int dir)
{
  Trace_Step step;
  step.pid = pid;
  step.pos = pos;
  step.dir = dir;
  chain->push_back(step);
}

/**/

static void append_chain(Trace_Chain * chain,const Trace_Chain & other,int shift)
{
  for (size_t i = 0; i < other.size();i++)
    push_step(chain,other[i].pid,other[i].pos+shift,other[i].dir);
}

/**/

bool Proof_Trace::cancel(Trace_Word * lhs,Trace_Word * rhs,Trace_Chain * chain,
                         size_t leading,size_t trailing,std::string * why) const
{
  /* lhs = P*a*S and rhs = P*b*S, and chain derives one from the other.
     Make it a derivation of a = b instead */
  if (leading + trailing > lhs->size() || leading + trailing > rhs->size())
  {
    *why = "cancellation longer than the words";
    return false;
  }
  if (trailing)
  {
    /* s(1)...s(k) = S and s(i)*w(i) = 1. Insert s(1)*w(1), then s(2)*w(2)
       inside it and so on, run the chain, and take them out again */
    size_t x = lhs->size() - trailing;
    size_t y = rhs->size() - trailing;
    Trace_Chain answer;
    std::vector<int> pid(trailing);
    for (size_t i = 0; i < trailing;i++)
    {
      Ordinal s = (*lhs)[x+i];
      if ((*rhs)[y+i] != s)
      {
        *why = "cancelled suffix differs";
        return false;
      }
      pid[i] = right_inverse[s];
      if (!pid[i])
      {
        *why = "no right inverse equation for a cancelled generator";
        return false;
      }
      push_step(&answer,pid[i],int(x+i),-1);
    }
    append_chain(&answer,*chain,0);
    for (size_t i = trailing; i-- > 0;)
      push_step(&answer,pid[i],int(y+i),1);
    chain->swap(answer);
    lhs->resize(x);
    rhs->resize(y);
  }
  if (leading)
  {
    /* p(1)...p(k) = P and w(i)*p(i) = 1. Insert w(k)*p(k) in front, then
       w(k-1)*p(k-1) between w(k) and p(k) and so on, run the chain shifted
       along, and take them out again */
    std::vector<int> pid(leading);
    std::vector<size_t> at(leading+1);
    at[leading] = 0;
    for (size_t i = leading; i-- > 0;)
    {
      Ordinal p = (*lhs)[i];
      if ((*rhs)[i] != p)
      {
        *why = "cancelled prefix differs";
        return false;
      }
      pid[i] = left_inverse[p];
      if (!pid[i])
      {
        *why = "no left inverse equation for a cancelled generator";
        return false;
      }
      /* at[i] is where w(i)*p(i) goes: after w(k)...w(i+1) */
      at[i] = at[i+1] + (i+1 < leading ? lhs_of[pid[i+1]].size()-1 : 0);
    }
    Trace_Chain answer;
    for (size_t i = leading; i-- > 0;)
      push_step(&answer,pid[i],int(at[i]),-1);
    size_t w_length = at[0] + lhs_of[pid[0]].size()-1;
    append_chain(&answer,*chain,int(w_length));
    for (size_t i = 0; i < leading;i++)
      push_step(&answer,pid[i],int(at[i]),1);
    chain->swap(answer);
    lhs->erase(lhs->begin(),lhs->begin()+leading);
    rhs->erase(rhs->begin(),rhs->begin()+leading);
  }
  return true;
}

/**/

void Proof_Trace::complain(const char * what,const std::string & why) const
{
  if (debug)
    fprintf(stderr,"#trace: %s: %s\n",what,why.c_str());
}

/**/

int Proof_Trace::axiom(const Word & lhs,const Word & rhs)
{
  Trace_Word l,r;
  read_word(&l,lhs);
  read_word(&r,rhs);
  return write("axiom",l,r,Trace_Chain());
}

/**/

int Proof_Trace::emit(const char * origin,const Word & lhs,const Word & rhs,
                      const Trace_Chain & chain)
{
  Trace_Word l,r;
  read_word(&l,lhs);
  read_word(&r,rhs);
  std::string why;
  if (!replay(l,chain,r,&why))
  {
    complain(origin,why);
    Trace_Chain hole;
    push_step(&hole,0,0,1);
    return write(origin,l,r,hole);
  }
  return write(origin,l,r,chain);
}

/**/

static const char * derivation_name(Basic_Derivation_Type type)
{
  switch (type)
  {
    case BDT_Trivial: return "triv";
    case BDT_Axiom: return "given";
    case BDT_Known: return "known";
    case BDT_Pool: return "pool";
    case BDT_Right_Conjugate: return "rconj";
    case BDT_Left_Conjugate: return "lconj";
    case BDT_Conjunction: return "ovl";
    case BDT_Equal_RHS: return "eqrhs";
    case BDT_Reverse: return "reverse";
    case BDT_Update: return "update";
    case BDT_Half_Difference_1: return "halfdiff";
    case BDT_Half_Difference_2: return "halfdiff";
    case BDT_Inversion_Closure: return "invclosure";
    case BDT_Diff_Reduction: return "wdiff";
    case BDT_Double: return "double";
    case BDT_Equal_LHS: return "eqlhs";
    case BDT_Relator_Multiply: return "relmult";
    case BDT_Normal_Closure: return "nclosure";
    case BDT_Pool_Conjugate: return "pconj";
    case BDT_Unspecified:
    default: return "unknown";
  }
}

/**/

void Proof_Trace::begin(Working_Equation & we)
{
  /* Called the first time the equation is logged, when its words are the
     ones the derivation type talks about. Works out the chain from the
     left hand side to the right hand side at that point. */
  Full_Derivation & d = we.derivation;
  const Node_Manager & nm = we.nm;
  Trace_Chain & base = d.trace_base;
  char origin[64];
  base.clear();
  d.trace_ops.clear();
  d.traced = true;
  strcpy(origin,derivation_name(d.type));

  switch (d.type)
  {
    case BDT_Axiom:
      {
        /* An axiom of the presentation gets its own line; what is stored
           is that seen through normalisation. MAF also re-adds its
           "polished" axioms, reduced forms of the axioms it derived, as
           axioms when it restarts: those must be cited, not assumed */
        Trace_Word l,r;
        read_word(&l,d.original_lhs);
        read_word(&r,d.original_rhs);
        int ax = content(l,r);
        if (!ax && given.count(key(l,r)))
          ax = axiom(d.original_lhs,d.original_rhs);
        push_step(&base,ax,0,1);
        sprintf(origin,"norm(%d)",ax);
      }
      break;
    case BDT_Known:
      if (!d.data.overlap.e1)
      {
        /* an equation read back in by -resume: all we can do is take it as
           given */
        resumed = true;
        int ax = axiom(d.original_lhs,d.original_rhs);
        push_step(&base,ax,0,1);
        sprintf(origin,"norm(%d)",ax);
      }
      else
      {
        /* lhs(e) to its right hand side, which is the raw one with -v -ve
           and otherwise the fully reduced one: as far along the chain of
           right hand sides as the words say */
        node_chain(nm,d.data.overlap.e1,0,&base);
        Trace_Word l,r;
        read_word(&l,d.original_lhs);
        read_word(&r,d.original_rhs);
        std::string why;
        for (size_t k = 1; k < base.size();k++)
          if (replay(l,Trace_Chain(base.begin(),base.begin()+k),r,&why))
          {
            base.resize(k);
            break;
          }
        sprintf(origin,"simp(%d)",base[0].pid);
      }
      break;
    case BDT_Trivial:
      break;
    case BDT_Conjunction:
      {
        /* the overlap of lhs(e1) and lhs(e2) rewritten both ways: back
           from rhs(e1)*tail to the overlap word, then on to head*rhs(e2) */
        node_chain(nm,d.data.overlap.e1,0,&base);
        reverse_chain(&base);
        int p1 = base.empty() ? 0 : base.back().pid;
        size_t n1 = base.size();
        node_chain(nm,d.data.overlap.e2,d.data.overlap.offset,&base);
        sprintf(origin,"ovl(%d,%d)",p1,base[n1].pid);
      }
      break;
    case BDT_Pool:
      {
        /* the pool equation with the words we were given */
        Trace_Word l,r;
        read_word(&l,d.original_lhs);
        read_word(&r,d.original_rhs);
        int p = content(l,r);
        push_step(&base,p,0,1);
        sprintf(origin,"simp(%d)",p);
      }
      break;
    case BDT_Equal_LHS:
      /* two right hand sides of the same left hand side: back from the old
         one, then forward to the new one */
      push_step(&base,d.tpid[0],0,-1);
      push_step(&base,d.tpid[1],0,1);
      sprintf(origin,"eqlhs(%d,%d)",d.tpid[0],d.tpid[1]);
      break;
    case BDT_Right_Conjugate:
    case BDT_Left_Conjugate:
    case BDT_Equal_RHS:
    case BDT_Pool_Conjugate:
      {
        Trace_Word l,r;
        read_word(&l,d.original_lhs);
        read_word(&r,d.original_rhs);
        std::vector<Trace_Chain> candidates;
        int p1 = 0,p2 = 0;
        conjugates(nm,d,l,r,&candidates,&p1,&p2);
        if (!justify(l,r,candidates,&base))
          push_step(&base,0,0,1);
        if (p2)
          sprintf(origin,"%s(%d,%d)",derivation_name(d.type),p1,p2);
        else
          sprintf(origin,"%s(%d)",derivation_name(d.type),p1);
      }
      break;
    default:
      /* Something we have no way to turn into rewriting steps, such as an
         equation proved by word-differences */
      push_step(&base,0,0,1);
      break;
  }
  d.trace_origin = origin;
}

/**/

bool Proof_Trace::justify(const Trace_Word & lhs,const Trace_Word & rhs,
                          const std::vector<Trace_Chain> & candidates,
                          Trace_Chain * chain) const
{
  /* The first candidate that derives rhs from lhs, either way round */
  std::string why;
  for (size_t i = 0; i < candidates.size();i++)
  {
    if (replay(lhs,candidates[i],rhs,&why))
    {
      *chain = candidates[i];
      return true;
    }
    Trace_Chain reversed = candidates[i];
    reverse_chain(&reversed);
    if (replay(lhs,reversed,rhs,&why))
    {
      *chain = reversed;
      return true;
    }
  }
  return false;
}

/**/

void Proof_Trace::conjugates(const Node_Manager & nm,const Full_Derivation & d,
                             const Trace_Word & l,const Trace_Word & r,
                             std::vector<Trace_Chain> * candidates,
                             int * p1,int * p2)
{
  /* Candidate derivations for the equations MAF deduces by multiplying
     an equation by a generator and cancelling, which need an equation
     g*w = 1 or w*g = 1 besides the ones the derivation names. Which way
     round the result is, and whether the right hand side has been reduced
     further than the cited equation's, varies between the places these
     come from, so justify() tries each. */
  Trace_Chain c;
  switch (d.type)
  {
    case BDT_Right_Conjugate:
      {
        /* e: a*x = b gives b*X = a: b*X -> a*x*X -> a */
        int p = d.tpid[0] ? d.tpid[0] : node_equation(nm,d.data.conjugate.e1);
        *p1 = p;
        if (!p || lhs_of[p].empty())
          return;
        const Trace_Word & ax = lhs_of[p];
        Trace_Word xX(2);
        xX[0] = ax.back();
        xX[1] = inverse_of[ax.back()];
        int inv = xX[1] == INVALID_SYMBOL ? 0 : content(xX,Trace_Word());
        push_step(&c,p,0,-1);
        push_step(&c,inv,int(ax.size()-1),1);
        candidates->push_back(c);
      }
      break;
    case BDT_Left_Conjugate:
      {
        Ordinal g = d.data.conjugate.g;
        Node_ID e = d.data.conjugate.e1;
        int p = node_equation(nm,e);
        *p1 = p;
        if (!p || lhs_of[p].empty())
          return;
        Trace_Chain e_at_1;
        node_chain(nm,e,1,&e_at_1);
        /* g*u = g*v from u = v */
        push_step(&c,p,1,1);
        candidates->push_back(c);
        candidates->push_back(e_at_1);
        /* e: x*u = v with g = inv(x) gives u = g*v: u -> g*x*u -> g*v */
        Trace_Word gx(2);
        gx[0] = g;
        gx[1] = lhs_of[p][0];
        int inv = content(gx,Trace_Word());
        c.clear();
        push_step(&c,inv,0,-1);
        push_step(&c,p,1,1);
        candidates->push_back(c);
        c.resize(1);
        append_chain(&c,e_at_1,0);
        candidates->push_back(c);
      }
      break;
    case BDT_Equal_RHS:
      {
        /* Two equations with a common side, after rewriting, and perhaps a
           common last generator y that can be cancelled with y*w = 1:
           A -> A*y*w -> ... -> B*y*w -> B */
        int p = d.tpid[0] ? d.tpid[0] : node_equation(nm,d.data.overlap.e1);
        int q = d.tpid[1] ? d.tpid[1] : node_equation(nm,d.data.overlap.e2);
        *p1 = p;
        *p2 = q;
        if (!p || !q)
          return;
        std::vector<Trace_Chain> middles;
        int e[2] = {p,q};
        for (int i = 0; i < 2;i++)
          for (int di = -1; di <= 1;di += 2)
          {
            c.clear();
            push_step(&c,e[i],0,di);
            middles.push_back(c);
            for (int dj = -1; dj <= 1;dj += 2)
            {
              c.resize(1);
              push_step(&c,e[1-i],0,dj);
              middles.push_back(c);
            }
          }
        for (size_t i = 0; i < middles.size();i++)
          candidates->push_back(middles[i]);
        std::vector<Ordinal> ys;
        const Trace_Word * words[4] = {&lhs_of[p],&lhs_of[q],&rhs_of[p],&rhs_of[q]};
        for (int i = 0; i < 4;i++)
          if (!words[i]->empty() &&
              std::find(ys.begin(),ys.end(),words[i]->back()) == ys.end())
            ys.push_back(words[i]->back());
        for (size_t k = 0; k < ys.size();k++)
        {
          int w = right_inverse[ys[k]];
          if (!w)
            continue;
          for (size_t i = 0; i < middles.size();i++)
          {
            c.clear();
            push_step(&c,w,int(l.size()),-1);
            append_chain(&c,middles[i],0);
            push_step(&c,w,int(r.size()),1);
            candidates->push_back(c);
            c.clear();
            push_step(&c,w,int(r.size()),-1);
            append_chain(&c,middles[i],0);
            push_step(&c,w,int(l.size()),1);
            candidates->push_back(c);
          }
        }
      }
      break;
    case BDT_Pool_Conjugate:
      {
        /* pool equation u*x = v gives u = v*X: u -> u*x*X -> v*X. One of
           our words is u and the other v*X */
        for (int side = 0; side < 2;side++)
        {
          const Trace_Word & u = side ? r : l;
          const Trace_Word & vX = side ? l : r;
          if (vX.empty() || inverse_of[vX.back()] == INVALID_SYMBOL)
            continue;
          Ordinal x = inverse_of[vX.back()];
          Trace_Word ux(u);
          ux.push_back(x);
          Trace_Word v(vX.begin(),vX.end()-1);
          int p = content(ux,v);
          Trace_Word xX(2);
          xX[0] = x;
          xX[1] = vX.back();
          int inv = content(xX,Trace_Word());
          if (!p || !inv)
            continue;
          *p1 = p;
          c.clear();
          push_step(&c,inv,int(u.size()),-1);
          push_step(&c,p,0,1);
          if (side)
            reverse_chain(&c);
          candidates->push_back(c);
        }
      }
      break;
    default:
      break;
  }
}

/**/

int Proof_Trace::emit_equation(Working_Equation & we)
{
  Full_Derivation & d = we.derivation;
  Trace_Word lhs,rhs,a,b;
  read_word(&lhs,we.lhs_word());
  read_word(&rhs,we.rhs_word());
  std::string why;
  Trace_Chain chain;
  bool ok = d.traced && d.started;

  if (!ok)
    why = "nothing was recorded";
  else
  {
    read_word(&a,d.original_lhs);
    read_word(&b,d.original_rhs);
    chain = d.trace_base;
    ok = replay(a,chain,b,&why);
    if (!ok)
      why = "start: " + why;
  }

  size_t failed_op = 0;
  for (size_t i = 0; ok && i < d.trace_ops.size();i++)
  {
    const Trace_Op & op = d.trace_ops[i];
    failed_op = i;
    Trace_Step step;
    step.pid = op.pid;
    step.pos = op.pos;
    switch (op.type)
    {
      /* MAF's record of reductions can carry an entry over from a word
         that has since been replaced or moved to the other side. A rewrite
         that does not apply where it says is such a leftover, and is
         skipped; whatever is skipped, the replay still has to end at the
         equation MAF has, or the derivation is written as a hole */
      case DST_Left_Reduction:
        step.dir = 1;
        if (step.pid && !apply(&a,step,&why))
          break;
        step.dir = -1;
        chain.insert(chain.begin(),step);
        if (!step.pid)
        {
          ok = false;
          why = "a reduction by an equation the trace does not have";
        }
        break;
      case DST_Right_Reduction:
        step.dir = 1;
        if (step.pid && !apply(&b,step,&why))
          break;
        chain.push_back(step);
        if (!step.pid)
        {
          ok = false;
          why = "a reduction by an equation the trace does not have";
        }
        break;
      case DST_Swap:
        a.swap(b);
        reverse_chain(&chain);
        break;
      case DST_Cancellation:
        ok = cancel(&a,&b,&chain,op.leading,op.trailing,&why);
        break;
      case DST_Right_Balance:
        {
          /* a = t*x becomes t, b becomes b*g, g = inv(x): t -> t*x*g */
          Trace_Word xg(2);
          if (a.empty())
          {
            ok = false;
            why = "balancing an empty word";
            break;
          }
          xg[0] = a.back();
          xg[1] = op.g;
          int p = content(xg,Trace_Word());
          a.pop_back();
          Trace_Chain answer;
          push_step(&answer,p,int(a.size()),-1);
          append_chain(&answer,chain,0);
          chain.swap(answer);
          b.push_back(op.g);
          if (!p)
          {
            ok = false;
            why = "no equation x*inv(x) = 1 for right balancing";
          }
        }
        break;
      case DST_Left_Balance:
        {
          /* a = x*t becomes t, b becomes g*b, g = inv(x): t -> g*x*t */
          Trace_Word gx(2);
          if (a.empty())
          {
            ok = false;
            why = "balancing an empty word";
            break;
          }
          gx[0] = op.g;
          gx[1] = a.front();
          int p = content(gx,Trace_Word());
          a.erase(a.begin());
          Trace_Chain answer;
          push_step(&answer,p,0,-1);
          append_chain(&answer,chain,1);
          chain.swap(answer);
          b.insert(b.begin(),op.g);
          if (!p)
          {
            ok = false;
            why = "no equation inv(x)*x = 1 for left balancing";
          }
        }
        break;
      case DST_Unbalance:
        {
          /* b = c*y becomes c, a becomes a*g, g = inv(y): c*y*g -> c */
          Trace_Word yg(2);
          if (b.empty())
          {
            ok = false;
            why = "unbalancing an empty word";
            break;
          }
          yg[0] = b.back();
          yg[1] = op.g;
          int p = content(yg,Trace_Word());
          b.pop_back();
          a.push_back(op.g);
          push_step(&chain,p,int(b.size()),1);
          if (!p)
          {
            ok = false;
            why = "no equation y*inv(y) = 1 for unbalancing";
          }
        }
        break;
      case DST_Right_Free_Reduction:
        {
          /* inv(g)*g removed from the end of the rhs */
          Trace_Word hg(2);
          hg[0] = inverse_of[op.g];
          hg[1] = op.g;
          int p = content(hg,Trace_Word());
          if (b.size() < 2 || b[b.size()-2] != hg[0] || b.back() != hg[1] || !p)
          {
            ok = false;
            why = "free reduction not justified";
            break;
          }
          step.pid = p;
          step.pos = int(b.size()-2);
          step.dir = 1;
          b.resize(b.size()-2);
          chain.push_back(step);
        }
        break;
      default:
        ok = false;
        why = "an operation the trace cannot express";
        break;
    }
  }

  if (ok && (a != lhs || b != rhs))
  {
    ok = false;
    why = "the recorded operations do not produce the equation";
  }
  const char * origin = d.traced ? d.trace_origin.c_str() : "unrecorded";
  if (!ok)
  {
    complain(origin,why);
    if (debug)
    {
      std::string text;
      char buffer[64];
      for (size_t i = 0; i < d.trace_ops.size();i++)
      {
        sprintf(buffer,"%s%d(%d@%d g%d %d/%d)",i == failed_op ? " *" : " ",
                d.trace_ops[i].type,d.trace_ops[i].pid,d.trace_ops[i].pos,
                d.trace_ops[i].g,d.trace_ops[i].leading,d.trace_ops[i].trailing);
        text += buffer;
      }
      complain("  ops",text);
      Trace_Word l0,r0;
      read_word(&l0,d.original_lhs);
      read_word(&r0,d.original_rhs);
      const Trace_Word * w[6] = {&l0,&r0,&a,&b,&lhs,&rhs};
      const char * names[6] = {"  start lhs","  start rhs","  got lhs","  got rhs",
                               "  want lhs","  want rhs"};
      for (int i = 0; i < 6;i++)
      {
        text.clear();
        for (size_t j = 0; j < w[i]->size();j++)
        {
          sprintf(buffer," %d",(*w[i])[j]+1);
          text += buffer;
        }
        complain(names[i],text);
      }
    }
    chain.clear();
    push_step(&chain,0,0,1);
  }
  return write(origin,lhs,rhs,chain);
}

/**/

void Proof_Trace::final_system(const Rewriting_System & rws,const Node_Manager & nm,
                               bool lostinfo)
{
  /* The rules written out are each the lhs of an equation node with its
     fully reduced right hand side; where that is not the right hand side
     the node's equation has, derive it */
  final_pids.clear();
  Element_Count count = rws.equation_count();
  Ordinal_Word lhs_word(nm.maf.alphabet);
  Ordinal_Word rhs_word(nm.maf.alphabet);
  for (Element_ID i = 1; i < count;i++)
  {
    rws.read_lhs(&lhs_word,i);
    rws.read_rhs(&rhs_word,i);
    Trace_Word l,r;
    read_word(&l,lhs_word);
    read_word(&r,rhs_word);
    int pid = content(l,r);
    if (!pid)
    {
      Node_Reference n = nm.get_state(lhs_word.buffer(),lhs_word.length());
      Trace_Chain chain;
      char origin[64];
      if (!n.is_null() && n->is_equation() && n->length(nm) == lhs_word.length())
      {
        node_chain(nm,n,0,&chain);
        sprintf(origin,"simp(%d)",chain[0].pid);
      }
      else
        strcpy(origin,"final");
      pid = emit(origin,lhs_word,rhs_word,chain);
    }
    final_pids.push_back(pid);
  }
  have_final = true;
  final_confluent = rws.is_confluent();
  final_lostinfo = lostinfo;
  fflush(file);
}
