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


// $Log: automata.cpp $
// Revision 1.15  2010/06/17 10:08:36Z  Alun
// Minor changes to option handling
// Revision 1.14  2010/06/10 17:38:16Z  Alun
// Fixes to ensure clean compilation with g++
// Revision 1.13  2010/06/10 13:56:50Z  Alun
// All tabs removed again
// Revision 1.12  2010/06/10 09:00:34Z  Alun_Williams
// Many options added or changed, especially new style strategy option
// Proper help finally added
// Revision 1.11  2009/11/10 08:34:43Z  Alun
// Command line option handling changed, + several new options added
// Revision 1.10  2009/10/13 23:07:58Z  Alun
// Compatibility of FSA utilities with KBMAG improved.
// Compatibility programs produce fewer FSA KBMAG would not (if any)
// Revision 1.9  2009/09/12 18:47:09Z  Alun
// printf() style functions cleaned up considerably
// code now compiles cleanly in g++ both 32-bit and 64-bit
// all tabs removed from source files
// Revision 1.8  2008/10/21 13:03:26Z  Alun
// Re-check in of source from web site after catastrophic loss of correct state
// of RCS archive.
// Revision 1.9  2008/10/21 14:03:26Z  Alun
// Early November checkin
// Completely reworked pool and Difference_Tracker
// Ported to Darwin
// Revision 1.8  2008/10/08 00:10:24Z  Alun
// force_multiplier option added
// Revision 1.7  2008/09/25 09:35:54Z  Alun
// Final version built using Indexer.
// Revision 1.5  2007/11/15 22:58:09Z  Alun
//
#include <time.h>
#include <string.h>
#include <signal.h>
#include <stdlib.h>
#include "maf.h"
#include "fsa.h"
#include "container.h"
#include "maf_sub.h"
#include "maf_so.h"
#include "mafctype.h"
#include "maf_trace.h"

#ifndef PROGRAM_NAME
#define PROGRAM_NAME "automata"
#endif

#ifndef COSETS
#define COSETS false
#endif

#ifndef NO_WD
#define NO_WD false
#endif

#ifndef WD
#define WD false
#endif

#ifndef NO_KB
#define NO_KB 0
#endif

#ifndef GA_FLAGS
#define GA_FLAGS GA_ALL
#endif

#ifndef VALIDATE
#define VALIDATE false
#endif

#ifndef EXCLUDE
#define EXCLUDE 0
#endif

#ifndef IS_KBPROG
#define IS_KBPROG 0
#endif

static MAF * global_maf;

void signal_handler(int /*signal_nr*/)
{
  if (global_maf)
    global_maf->give_up();
}

int main(int argc,char ** argv)
{
  int i = 1;
  char * group_filename = 0;
  char * subgroup_suffix = 0;
  Container & container = *MAF::create_container();
  Standard_Options so(container,SO_FSA_FORMAT);
  bool bad_usage = false;
  FSA_Simple * fsa = 0;
  unsigned flags = CFI_DEFAULT|CFI_ALLOW_CREATE|CFI_CREATE_RM;
#define cprintf container.error_output
  bool no_validate = false;
  bool cosets = COSETS;
  MAF::Options options;
  int exit_code = 0;
  unsigned exclude = EXCLUDE;
#if defined(IS_KBPROG) && IS_KBPROG
  bool both = false;
#endif

  options.no_differences = NO_WD;
  options.differences = WD;
  options.no_kb = NO_KB;
  options.force_differences = NO_KB;
  options.validate = VALIDATE;
  options.write_success = VALIDATE;
#ifdef EMULATE
  options.emulate_kbmag = true;
#endif
  options.is_kbprog = IS_KBPROG;

  while (i < argc && !bad_usage)
  {
    if (argv[i][0] == '-')
    {
      String arg = argv[i];
      if (arg.is_equal("-balancing"))
      {
        so.parse_natural(&options.right_balancing,argv[i+1],7,arg);
        i += 2;
      }
#if IS_KBPROG
      else if (arg.is_equal("-both"))
      {
        both = true;
        i++;
      }
#endif
      else if (arg.is_equal("-check"))
      {
        fsa = FSA_Factory::create(argv[i+1],&container);
        i += 2;
      }
      else if (so.present(arg,"-check_inverses"))
      {
        options.check_inverses = true;
        i++;
      }
      else if (arg.is_equal("-collapse"))
      {
        so.parse_natural(&options.collapse,argv[i+1],2,arg);
        i += 2;
      }
      else if (arg.is_equal("-confluent"))
      {
        options.assume_confluent = true;
        i++;
      }
      else if (arg.is_equal("-conjugation"))
      {
        so.parse_natural(&options.conjugation,argv[i+1],4,arg);
        i += 2;
      }
      else if (so.present(arg,"-consider_secondary"))
      {
        options.consider_secondary = true;
        i++;
      }
      else if (arg.is_equal("-cos"))
      {
        cosets = true;
        if (!no_validate)
          options.validate = true;
        i++;
      }
      else if (so.present(arg,"-dense_rm"))
      {
        options.dense_rm = true;
        i++;
      }
      else if (
#ifdef KBMAG_FINITE_INDEX
               arg.is_equal("-f") ||
#endif
               so.present(arg,"-detect_finite_index"))
      {
        options.detect_finite_index = 1;
        i++;
      }
      else if (arg.is_equal("-eliminate"))
      {
        options.eliminate = true;
        i++;
      }
      else if (so.present(arg,"-expand_all"))
      {
        options.expand_all = true;
        i++;
      }
      else if (so.present(arg,"-extended_consider"))
      {
        options.extended_consider = true;
        i++;
      }
      else if (arg.is_equal("-fast"))
      {
        options.fast = true;
        i++;
      }
      else if (arg.is_equal("-filters"))
      {
        so.parse_natural(&options.filters,argv[i+1],7,arg);
        i += 2;
      }
      else if (so.present(arg,"-force_differences"))
      {
        options.force_differences = true;
        i++;
      }
      else if (so.present(arg,"-force_multiplier"))
      {
        options.force_multiplier = true;
        i++;
      }
#if defined(EMULATE) && !IS_KBPROG
      else if (arg.is_equal("-gpminkb"))
      {
        exclude &= ~(GA_L1_ACCEPTOR|GA_DIFF2C|GA_DIFF1C|GA_MINKB);
        i++;
      }
#endif
      else if (so.present(arg,"-ignore_h_length"))
      {
        options.ignore_h_length = true;
        i++;
      }
      else if (so.present(arg,"-left_balancing"))
      {
        so.parse_natural(&options.left_balancing,argv[i+1],3,arg);
        i += 2;
      }
      else if (arg.is_equal("-m") || arg.is_equal("-me") || so.present(arg,"-max_equations"))
      {
        so.parse_natural(&options.max_equations,argv[i+1],0,arg);
        i += 2;
      }
      else if (arg.is_equal("-mo") || so.present(arg,"-max_overlap"))
      {
        so.parse_total_length(&options.max_overlap_length,argv[i+1],0,arg);
        i += 2;
      }
      else if (so.present(arg,"-max_time"))
      {
        so.parse_natural(&options.max_time,argv[i+1],0,arg);
        i += 2;
      }
      else if (arg.is_equal("-mt") || so.present(arg,"-min_time"))
      {
        so.parse_natural(&options.min_time,argv[i+1],0,arg);
        i += 2;
      }
      else if (so.present(arg,"-no_early_repair"))
      {
        options.no_early_repair = true;
        i++;
      }
      else if (so.present(arg,"-no_h"))
      {
        options.no_h = true;
        i++;
      }
      else if (so.present(arg,"-no_hd"))
      {
        options.no_half_differences = true;
        i++;
      }
      else if (so.present(arg,"-no_kb"))
      {
        options.no_kb = true;
        i++;
      }
      else if (so.present(arg,"-no_pool"))
      {
        options.no_pool = true;
        i++;
      }
      else if (so.present(arg,"-no_pool_below"))
      {
        so.parse_natural(&options.no_pool_below,argv[i+1],MAX_WORD,arg);
        i += 2;
      }
      else if (so.present(arg,"-no_prune"))
      {
        options.no_prune = true;
        i++;
      }
      else if (so.present(arg,"-no_throttle"))
      {
        options.no_throttle = true;
        i++;
      }
      else if (so.present(arg,"-no_unbalance"))
      {
        options.unbalance = false;
        i++;
      }
      else if (so.present(arg,"-no_validate"))
      {
        no_validate = true;
        options.validate = false;
        i++;
      }
      else if (so.present(arg,"-no_weak_acceptor"))
      {
        options.no_weak_acceptor = true;
        i++;
      }
#if !NO_KB
      else if (so.present(arg,"-no_wd"))
      {
        options.no_differences = true;
        options.differences = false;
        i++;
      }
#endif
      else if (arg.is_equal("-p") || arg.is_equal("-subpres"))
      {
        exclude &= ~GA_SUBPRES;
        i++;
      }
      else if (so.present(arg,"-partial_reductions"))
      {
        so.parse_natural(&options.partial_reductions,argv[i+1],2,arg);
        i += 2;
      }
      else if (so.present(arg,"-pool_above"))
      {
        so.parse_total_length(&options.pool_above,argv[i+1],0,arg);
        i += 2;
      }
      else if (so.present(arg,"-probe_style"))
      {
        so.parse_natural(&options.probe_style,argv[i+1],3,arg);
        i += 2;
      }
      else if (so.present(arg,"-prove_finite_index"))
      {
        options.detect_finite_index = 2;
        i++;
      }
      else if (arg.is_equal("-repeat"))
      {
        so.parse_natural(&options.repeat,argv[i+1],2,arg);
        i += 2;
      }
      else if (arg.is_equal("-resume") || arg.is_equal("-r") || arg.is_equal("-ro"))
      {
        flags |= CFI_RESUME;
        i++;
      }
      else if (arg.is_equal("-schreier"))
      {
        options.use_schreier_generators = true;
        i++;
      }
      else if (arg.is_equal("-secondary"))
      {
        so.parse_natural(&options.secondary,argv[i+1],4,arg);
        i += 2;
      }
      else if (arg.is_equal("-slow"))
      {
        options.slow = true;
        i++;
      }
      else if (so.present(arg,"-special_overlaps"))
      {
        so.parse_natural(&options.special_overlaps,argv[i+1],5,arg);
        i += 2;
      }
      else if (arg.is_equal("-strategy"))
      {
        static const char * strategy_name[] =
        {
          "naive",
          "vanilla",
          "easy",
          "quick",
          "short",
          "balanced",
          "wreath",
          "sparse",
          "era",
          "long",
          0
        };
        bool found = false;
        String arg = argv[i+1];
        String_Length mod_pos = 0;
        bool modified = arg.find('/',&mod_pos);
        Letter * sname = modified ? arg.clone_substring(0,mod_pos) : arg.clone();
        if (*sname)
        {
          unsigned local;
          for (local = 0;strategy_name[local]!=0;local++)
            if (so.present(sname,strategy_name[local]))
            {
              found = true;
              break;
            }
          if (found)
            options.set_strategy(MAF::MAF_Strategy(local));
          else
          {
            cprintf("Bad strategy name %s specified for -strategy\n",sname);
            bad_usage = true;
          }
        }
        else
          options.strategy = MSF_RESPECT_RIGHT;
        delete sname;
        if (modified)
        {
          for (const Letter * s = arg + mod_pos+1;*s;s++)
          {
            Letter c = char_classification.to_upper_case(*s);
            switch (c)
            {
              case 'A':
                options.strategy |= MSF_AGGRESSIVE_DISCARD;
                break;
              case 'B':
                options.strategy |= MSF_LEFT_EXPAND;
                break;
              case 'C':
                options.strategy |= MSF_CLOSED;
                break;
              case 'D':
                options.strategy |= MSF_DEEP_RIGHT;
                break;
              case 'E':
                options.strategy |= MSF_USE_ERAS;
                break;
              case 'F':
                options.strategy &= ~MSF_RESPECT_RIGHT;
                break;
              case 'I':
                options.no_deductions = false;
                break;
              case 'J':
                options.strategy |= MSF_CONJUGATE_FIRST;
                break;
              case 'L':
                options.strategy |= MSF_NO_FAVOUR_SHORT;
                break;
              case 'M':
                options.strategy &= ~MSF_DEEP_RIGHT;
                break;
              case 'N':
                options.strategy &= ~MSF_NO_FAVOUR_SHORT;
                break;
              case 'O':
                options.strategy &= ~MSF_CLOSED;
                break;
              case 'P':
                options.no_pool = false;
                break;
              case 'S':
                options.strategy |= MSF_SELECTIVE_PROBE;
                break;
              case 'T':
                options.strategy &= ~MSF_USE_ERAS;
                break;
              case 'U':
                options.no_deductions = true;
                break;
             default:
               cprintf("Unrecognised option '%c' found in -strategy string\n"
                       "(parsing %s)\n"
                       "(at      %*s)\n",
                       c,argv[i+1],int(s - arg+1),"^");
               bad_usage = true;
               break;
            }
          }
        }
        i += 2;
      }
      else if (so.present(arg,"-swap_bad"))
      {
        options.swap_bad = true;
        i++;
      }
      else if (arg.is_equal("-tight"))
      {
        options.tight = true;
        i++;
      }
      else if (arg.is_equal("-timeout"))
      {
        so.parse_natural(&options.timeout,argv[i+1],0,arg);
        i += 2;
      }
      else if (arg.is_equal("-trace"))
      {
        if (i+1 >= argc)
          bad_usage = true;
        else
          Proof_Trace::set_filename(argv[i+1]);
        i += 2;
      }
      else if (arg.is_equal("-validate"))
      {
        options.validate = true;
        i++;
      }
      else if (arg.is_equal("-ve"))
      {
        options.log_flags |= LOG_EQUATIONS;
        i++;
      }
      else if (arg.is_equal("-vwd"))
      {
        options.log_flags |= LOG_WORD_DIFFERENCES;
        i++;
      }
      else if (arg.is_equal("-wd"))
      {
        options.differences = true;
        options.no_differences = false;
        i++;
      }
      else if (so.present(arg,"-weed_secondary"))
      {
        options.weed_secondary = true;
        i++;
      }
      else if (so.present(arg,"-work_order"))
      {
        so.parse_natural(&options.expansion_order,argv[i+1],6,arg);
        i += 2;
      }
      else if (arg.is_equal("-l") ||
               arg.is_equal("-large") ||
               arg.is_equal("-d") ||
               arg.is_equal("-d1") ||
               arg.is_equal("-diff1") ||
               arg.is_equal("-f"))
        i++;  // for KBMAG compatibility
      else if (arg.is_equal("-cn") ||
               arg.is_equal("-hf") ||
               arg.is_equal("-mrl") ||
               arg.is_equal("-ms") ||
               arg.is_equal("-mwd") ||
               arg.is_equal("-sort") ||
               arg.is_equal("-t"))
      {
        i += 2;
        if (i > argc)
          bad_usage = true;
      }
      else if (arg.is_equal("-rk"))
      {
        i += 3;
        if (i > argc)
          bad_usage = true;
      }
      else if (!so.recognised(argv,i))
      {
        bad_usage = true;
        cprintf("Unrecognised command line option: %s\n",argv[i]);
      }
    }
    else
    {
      if (group_filename == 0)
        group_filename = argv[i];
      else if (cosets && subgroup_suffix == 0)
        subgroup_suffix = argv[i];
      else
        bad_usage = true;
      i++;
    }
  }
  if (options.log_flags & LOG_EQUATIONS && so.verbose)
    options.log_flags |= LOG_DERIVATIONS;
  if (Proof_Trace::wanted() && !bad_usage)
  {
    if (cosets)
    {
      cprintf("Proof tracing is not supported for coset systems.\n");
      return 1;
    }
    proof_trace = new Proof_Trace(Proof_Trace::wanted());
    if (flags & CFI_RESUME)
      proof_trace->set_resumed();
    options.log_flags |= LOG_TRACE;
  }
  container.set_log_level(so.log_level);
  if (group_filename && !bad_usage)
  {
    time_t now = time(0);
    MAF * maf = 0;
    options.log_level = so.log_level;

    options.log_level = so.log_level;
    options.tight = options.tight || options.assume_confluent || options.no_differences;
    options.assume_confluent |= options.no_differences;
#if IS_KBPROG
    if (!options.no_differences && !both)
      exclude |= GA_RWS;
#endif
    maf = MAF::create_from_input(cosets,group_filename,subgroup_suffix,
                                 &container,flags,&options);
    if (fsa)
      maf->set_validation_fsa(fsa);
    signal(SIGINT,signal_handler);
    signal(SIGTERM,signal_handler);
    global_maf = maf;
#ifdef GPMAKEFSA
    if (maf->load_fsas(GA_DIFF2))
    {
      maf->import_difference_machine(*maf->fsas.diff2);
      maf->grow_automata(0,GA_FLAGS,0,exclude);
    }
    else
    {
      container.error_output("Unable to load required word-difference"
                             " machine!\n");
      exit_code = 1;
    }
#else
    maf->grow_automata(0,GA_FLAGS,0,exclude);
#endif
    if (maf->aborting)
      exit_code = 2;
    if (proof_trace)
    {
      delete proof_trace; // writes the list of final equations
      proof_trace = 0;
    }
    container.progress(1,"Elapsed time %ld\n",long(time(0) - now));
    global_maf = 0;
    delete maf;
  }
  else
  {
    cprintf("Usage:" PROGRAM_NAME " [options] rwsname [-cos [subsuffix]]\n"
            "where rwsname is a GASP rewriting system, and, if the -cos"
            " option is used,\nfilename.subsuffix is substructure file.\n"
            "The most important command line options are summarised below,"
            " but for full\ndetails refer to the MAF documentation.\n"
#if !defined(EMULATE) && !defined(GPMAKEFSA)
"-nowd                Only compute automata not requiring word-differences to\n"
"                     be computed. This option is useful if you are only\n"
"                     interested in finding a confluent rewriting system.\n"
"-confluent           Look initially only for a confluent rewriting system, but\n"
"                     if one is found also construct automatic structures if\n"
"                     possible.\n"
"-nokb                Attempt to compute an automatic structure without using\n"
"                     Knuth-Bendix procedure first. This is is best for many\n"
"                     automatic groups.\n"
"-force_differences   Favour search for word-differences even if there are very\n"
"                     many of them. This option is useful if you want to do a\n"
"                     preliminary run to decide whether or not to look for an\n"
"                     automatic structure.\n"
#endif

"-trace file          Write a proof trace to file: each equation found, with the\n"
"                     rewriting steps that derive it, in the format of KBMAG's\n"
"                     kbprog -trace. Not for coset systems.\n"
#if !defined(EMULATE)
"-validate            Enable or disable axiom-checking for automatic structures.\n"
"-no_validate         The default behaviour is to check axioms only for coset\n"
"                     systems.\n"
#endif

#ifndef GPMAKEFSA
"-strategy name       Selects a predefined strategy. Usually MAF will select a\n"
"                     good strategy automatically, but if the input file has\n"
"                     a mix of long and short equations the \"long\" or \"era\"\n"
"                     strategies may work better. For simple input files \"easy\" or\n"
"                     \"quick\" strategies may sometimes work more quickly than the\n"
"                     default.\n\n"
"For coset systems:\n"
"\n"
"-prove_finite_index  Attempt to prove the subgroup has finite index, and if\n"
"                     successful output the rewriting system as soon as words\n"
"                     can be reduced to the minimal coset representative. The\n"
"                     -detect_finite_index option has a similar effect, but does\n"
"                     not modify the strategy.\n"
"-no_h                Ignore the H-equations when performing Knuth-Bendix.\n"
"-ignore_h_length     Ignore the H-word part of mixed words when deciding\n"
"                     whether to pool long equations.\n"
#endif

            );
    exit_code = 1;
  }

  delete &container;
  return exit_code;
}

