#Copyright 2008,2009,2010 Alun Williams
#This file is part of MAF.
#MAF is free software: you can redistribute it and/or modify it
#under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#MAF is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#You should have received a copy of the GNU General Public License
#along with MAF.  If not, see <http://www.gnu.org/licenses/>.

#$Log: mafu64.mak $
#Revision 1.6  2011/06/11 11:54:29Z  Alun
#Added fsaseparate and gpxlatwa
#Revision 1.5  2010/06/10 18:07:23Z  Alun
#Fixes to ensure clean compilation with g++

#If any of the compilation flags cause trouble, then you can probably remove them
#If you are not using GCC compilers then make sure you enable whatever flag
#ensures char is unsigned by default.
CPPFLAGS = -fPIC -m64 -Os -Wall -Wextra -Wno-switch -Wno-reorder -Wno-char-subscripts  -funsigned-char -fpermissive -fno-rtti -fno-threadsafe-statics -fno-exceptions 
LINKER = g++ -m64 
BIN = ./../bin
O=o
LINK_EXTRA =

all: \
  $(BIN)/example \
  $(LIBS) \
  $(BIN)/automata \
  $(BIN)/autcos \
  $(BIN)/autgroup \
  $(BIN)/kbprog \
  $(BIN)/kbprogcos \
  $(BIN)/fsaand \
  $(BIN)/fsaandnot \
  $(BIN)/fsabfs \
  $(BIN)/fsacartesian \
  $(BIN)/fsacompose \
  $(BIN)/fsaconcat \
  $(BIN)/fsacount \
  $(BIN)/fsacut \
  $(BIN)/fsadiagonal \
  $(BIN)/fsaenumerate \
  $(BIN)/fsaexists \
  $(BIN)/fsafl \
  $(BIN)/fsakernel \
  $(BIN)/fsalequal \
  $(BIN)/fsamerge \
  $(BIN)/fsamin \
  $(BIN)/fsan \
  $(BIN)/fsanot \
  $(BIN)/fsaor \
  $(BIN)/fsaprint \
  $(BIN)/fsaproduct \
  $(BIN)/fsaprune \
  $(BIN)/fsaread \
  $(BIN)/fsareverse \
  $(BIN)/fsapad \
  $(BIN)/fsaseparate \
  $(BIN)/fsashortlex \
  $(BIN)/fsastar \
  $(BIN)/fsaswapcoords \
  $(BIN)/fsatrim \
  $(BIN)/gpaxioms \
  $(BIN)/gpcclass \
  $(BIN)/gpcosets \
  $(BIN)/gpdifflabs \
  $(BIN)/gpgenmult2 \
  $(BIN)/gpgeowa \
  $(BIN)/gpmakefsa \
  $(BIN)/gpmigmdet \
  $(BIN)/gpmigenmult2 \
  $(BIN)/gpmimult \
  $(BIN)/gpmimult2 \
  $(BIN)/gpminkb \
  $(BIN)/gpmorphism \
  $(BIN)/gpmult \
  $(BIN)/gpmult2 \
  $(BIN)/gporder \
  $(BIN)/gpovlwa \
  $(BIN)/gpstabiliser \
  $(BIN)/gpsublowindex \
  $(BIN)/gpsubmake \
  $(BIN)/gpsubpres \
  $(BIN)/gpsubwa \
  $(BIN)/gptcenum \
  $(BIN)/gpvital \
  $(BIN)/gpxlatwa \
  $(BIN)/gpwa \
  $(BIN)/midfadeterminize \
  $(BIN)/reduce \
  $(BIN)/rwsprint \
  $(BIN)/makecosfile  \
  $(BIN)/isconjugate \
  $(BIN)/isnormal \
  $(BIN)/simplify

MAFLIB = \
  mafbase.$O \
  container.$O \
  platform.$O \
  stream.$O \
  maf.$O \
  present.$O \
  equation.$O \
  mafload.$O \
  mafread.$O \
  alphabet.$O \
  maf_avl.$O \
  maf_cfp.$O \
  maf_dr.$O \
  maf_dt.$O \
  maf_el.$O \
  maf_em.$O \
  maf_ew.$O \
  maf_jm.$O \
  maf_mult.$O \
  maf_nm.$O \
  maf_rm.$O \
  maf_rws.$O \
  maf_spl.$O \
  maf_so.$O \
  maf_ss.$O \
  maf_sub.$O \
  maf_subwa.$O \
  maf_we.$O \
  mafqueue.$O \
  mafctype.$O \
  mafword.$O \
  hash.$O \
  keyedfsa.$O \
  bitarray.$O \
  heap.$O \
  rubik.$O \
  fsa.$O \
  mafauto.$O  \
  mafconj.$O \
  mafcoset.$O \
  mafgeowa.$O \
  mafminkb.$O \
  nodelist.$O \
  mafnode.$O \
  ltfsa.$O \
  maf_tc.$O \
  relators.$O \
  lowindex.$O \
  subpres.$O \
  tietze.$O \
  x_to_str.$O

#If your system does not support dynamic libraries, or if you have trouble 
#getting MAF to work, then try linking against the static MAF library instead

#LIBS = maf.a
LIBS = $(BIN)/libmaf.dylib

maf.a : $(MAFLIB)
	rm -f $@
	ar rcs $@ $(MAFLIB)

#"make static" builds every program as a standalone executable in
#$(STATIC_BIN): linked statically against maf.a and the C and C++ runtime
#libraries, so that it needs neither libmaf.dylib nor anything else at run
#time and can be copied anywhere.
STATIC_BIN = ./../bin/static

static :
	-mkdir -p $(STATIC_BIN)
	$(MAKE) BIN=$(STATIC_BIN) LIBS=maf.a LINK_EXTRA=-static all

.PHONY : all static

$(BIN)/libmaf.dylib : $(MAFLIB)
	-mkdir $(BIN)
	g++ -m64 -shared -dynamiclib $(MAFLIB) -o $@

EXAMPLE = \
  example.$O \
  $(LIBS)

$(BIN)/example: $(EXAMPLE) 
	$(LINKER) -o $@ $(EXAMPLE) $(LINK_EXTRA)   

FSAREAD = \
  fsaread.$O \
  $(LIBS)

$(BIN)/fsaread: $(FSAREAD) 
	$(LINKER) -o $@ $(FSAREAD) $(LINK_EXTRA)   

REDUCE = \
  reduce.$O \
  $(LIBS)

$(BIN)/reduce: $(REDUCE) 
	$(LINKER) -o $@ $(REDUCE) $(LINK_EXTRA)   

AUTCOS = \
  autcos.$O \
  $(LIBS)  

$(BIN)/autcos: $(AUTCOS) 
	$(LINKER) -o $@ $(AUTCOS) $(LINK_EXTRA)   

AUTGROUP = \
  autgroup.$O \
  $(LIBS)

$(BIN)/autgroup: $(AUTGROUP) 
	$(LINKER) -o $@ $(AUTGROUP) $(LINK_EXTRA)   

AUTOMATA = \
  automata.$O \
  $(LIBS)

$(BIN)/automata: $(AUTOMATA) 
	$(LINKER) -o $@ $(AUTOMATA) $(LINK_EXTRA)   

KBPROG =  \
  kbprog.$O \
  $(LIBS)

$(BIN)/kbprog: $(KBPROG) 
	$(LINKER) -o $@ $(KBPROG) $(LINK_EXTRA)   

KBPROGCOS =  \
  kbprogcos.$O \
  $(LIBS)

$(BIN)/kbprogcos: $(KBPROGCOS) 
	$(LINKER) -o $@ $(KBPROGCOS) $(LINK_EXTRA)   

FSACONCAT = \
  fsaconcat.$O \
  $(LIBS)

$(BIN)/fsaconcat: $(FSACONCAT) 
	$(LINKER) -o $@ $(FSACONCAT) $(LINK_EXTRA)   

FSAPRODUCT = \
  fsaproduct.$O \
  $(LIBS)

$(BIN)/fsaproduct: $(FSAPRODUCT) 
	$(LINKER) -o $@ $(FSAPRODUCT) $(LINK_EXTRA)   

FSACOUNT = \
  fsacount.$O \
  $(LIBS)

$(BIN)/fsacount: $(FSACOUNT) 
	$(LINKER) -o $@ $(FSACOUNT) $(LINK_EXTRA)   

FSACUT = \
  fsacut.$O \
  $(LIBS)

$(BIN)/fsacut: $(FSACUT) 
	$(LINKER) -o $@ $(FSACUT) $(LINK_EXTRA)   

FSAPRINT = \
  fsaprint.$O \
  $(LIBS)

$(BIN)/fsaprint: $(FSAPRINT) 
	$(LINKER) -o $@ $(FSAPRINT) $(LINK_EXTRA)   

FSAAND = \
  fsaand.$O \
  $(LIBS)

$(BIN)/fsaand: $(FSAAND) 
	$(LINKER) -o $@ $(FSAAND) $(LINK_EXTRA)   

FSAOR = \
  fsaor.$O \
  $(LIBS)

$(BIN)/fsaor: $(FSAOR) 
	$(LINKER) -o $@ $(FSAOR) $(LINK_EXTRA)   

FSAANDNOT = \
  fsaandnot.$O \
  $(LIBS)

$(BIN)/fsaandnot: $(FSAANDNOT) 
	$(LINKER) -o $@ $(FSAANDNOT) $(LINK_EXTRA)   

FSABFS = \
  fsabfs.$O \
  $(LIBS)

$(BIN)/fsabfs: $(FSABFS) 
	$(LINKER) -o $@ $(FSABFS) $(LINK_EXTRA)   

FSADIAGONAL = \
  fsadiagonal.$O \
  $(LIBS)

$(BIN)/fsadiagonal: $(FSADIAGONAL) 
	$(LINKER) -o $@ $(FSADIAGONAL) $(LINK_EXTRA)   


FSAMERGE = \
  fsamerge.$O \
  $(LIBS)

$(BIN)/fsamerge: $(FSAMERGE) 
	$(LINKER) -o $@ $(FSAMERGE) $(LINK_EXTRA)   

FSAMIN = \
  fsamin.$O \
  $(LIBS)

$(BIN)/fsamin: $(FSAMIN) 
	$(LINKER) -o $@ $(FSAMIN) $(LINK_EXTRA)   

FSAPRUNE = \
  fsaprune.$O \
  $(LIBS)

$(BIN)/fsaprune: $(FSAPRUNE) 
	$(LINKER) -o $@ $(FSAPRUNE) $(LINK_EXTRA)   

FSAKERNEL = \
  fsakernel.$O \
  $(LIBS)

$(BIN)/fsakernel: $(FSAKERNEL) 
	$(LINKER) -o $@ $(FSAKERNEL) $(LINK_EXTRA)   

FSAENUMERATE = \
  fsaenumerate.$O \
  $(LIBS)

$(BIN)/fsaenumerate: $(FSAENUMERATE) 
	$(LINKER) -o $@ $(FSAENUMERATE) $(LINK_EXTRA)   

FSAEXISTS = \
  fsaexists.$O \
  $(LIBS)

$(BIN)/fsaexists: $(FSAEXISTS) 
	$(LINKER) -o $@ $(FSAEXISTS) $(LINK_EXTRA)   

FSAFL = \
  fsafl.$O \
  $(LIBS)

$(BIN)/fsafl: $(FSAFL) 
	$(LINKER) -o $@ $(FSAFL) $(LINK_EXTRA) 

FSAN = \
  fsan.$O \
  $(LIBS)

$(BIN)/fsan: $(FSAN) 
	$(LINKER) -o $@ $(FSAN) $(LINK_EXTRA)   

FSAPAD = \
  fsapad.$O \
  $(LIBS)

$(BIN)/fsapad: $(FSAPAD) 
	$(LINKER) -o $@ $(FSAPAD) $(LINK_EXTRA)   

FSASEPARATE = \
  fsaseparate.$O \
  $(LIBS)

$(BIN)/fsaseparate: $(FSASEPARATE) 
	$(LINKER) -o $@ $(FSASEPARATE) $(LINK_EXTRA)

FSASHORTLEX = \
  fsashortlex.$O \
  $(LIBS)

$(BIN)/fsashortlex: $(FSASHORTLEX) 
	$(LINKER) -o $@ $(FSASHORTLEX) $(LINK_EXTRA)   

FSASTAR = \
  fsastar.$O \
  $(LIBS)

$(BIN)/fsastar: $(FSASTAR) 
	$(LINKER) -o $@ $(FSASTAR) $(LINK_EXTRA)   

FSANOT = \
  fsanot.$O \
  $(LIBS)

$(BIN)/fsanot: $(FSANOT) 
	$(LINKER) -o $@ $(FSANOT) $(LINK_EXTRA)   

FSAREVERSE = \
  fsareverse.$O \
  $(LIBS)

$(BIN)/fsareverse: $(FSAREVERSE) 
	$(LINKER) -o $@ $(FSAREVERSE) $(LINK_EXTRA)   

FSASWAPCOORDS = \
  fsaswapcoords.$O \
  $(LIBS)

$(BIN)/fsaswapcoords: $(FSASWAPCOORDS) 
	$(LINKER) -o $@ $(FSASWAPCOORDS) $(LINK_EXTRA)   

FSATRIM = \
  fsatrim.$O \
  $(LIBS)

$(BIN)/fsatrim: $(FSATRIM) 
	$(LINKER) -o $@ $(FSATRIM) $(LINK_EXTRA)   

MIDFADETERMINIZE = \
  midfadeterminize.$O \
  $(LIBS)

$(BIN)/midfadeterminize: $(MIDFADETERMINIZE)
	$(LINKER) -o $@ $(MIDFADETERMINIZE) $(LINK_EXTRA)   

RWSPRINT = \
  rwsprint.$O \
  $(LIBS)
  
$(BIN)/rwsprint: $(RWSPRINT) 
	$(LINKER) -o $@ $(RWSPRINT) $(LINK_EXTRA)   

MAKECOSFILE = \
  makecos.$O \
  $(LIBS)

$(BIN)/makecosfile: $(MAKECOSFILE) 
	$(LINKER) -o $@ $(MAKECOSFILE) $(LINK_EXTRA)   

GPAXIOMS = \
  gpaxioms.$O \
  $(LIBS)

$(BIN)/gpaxioms: $(GPAXIOMS) 
	$(LINKER) -o $@ $(GPAXIOMS) $(LINK_EXTRA)   

GPCCLASS = \
  gpcclass.$O \
  $(LIBS)

$(BIN)/gpcclass: $(GPCCLASS) 
	$(LINKER) -o $@ $(GPCCLASS) $(LINK_EXTRA)   

GPCOSETS = \
  gpcosets.$O \
  $(LIBS)

$(BIN)/gpcosets: $(GPCOSETS) 
	$(LINKER) -o $@ $(GPCOSETS) $(LINK_EXTRA)   

GPMULT = \
  gpmult.$O \
  $(LIBS)

$(BIN)/gpmult: $(GPMULT) 
	$(LINKER) -o $@ $(GPMULT) $(LINK_EXTRA)   

GPMULT2 = \
  gpmult2.$O \
  $(LIBS)

$(BIN)/gpmult2: $(GPMULT2) 
	$(LINKER) -o $@ $(GPMULT2) $(LINK_EXTRA)   

GPMIMULT = \
  gpmimult.$O \
  $(LIBS)

$(BIN)/gpmimult: $(GPMIMULT) 
	$(LINKER) -o $@ $(GPMIMULT) $(LINK_EXTRA)   

GPMIMULT2 = \
  gpmimult2.$O \
  $(LIBS)

$(BIN)/gpmimult2: $(GPMIMULT2) 
	$(LINKER) -o $@ $(GPMIMULT2) $(LINK_EXTRA)   

GPMINKB = \
  gpminkb.$O \
  $(LIBS)

$(BIN)/gpminkb: $(GPMINKB) 
	$(LINKER) -o $@ $(GPMINKB) $(LINK_EXTRA)   

GPMIGMDET = \
  gpmigmdet.$O \
  $(LIBS)

$(BIN)/gpmigmdet: $(GPMIGMDET) 
	$(LINKER) -o $@ $(GPMIGMDET) $(LINK_EXTRA)   

GPMAKEFSA = \
  gpmakefsa.$O \
  $(LIBS)

$(BIN)/gpmakefsa: $(GPMAKEFSA) 
	$(LINKER) -o $@ $(GPMAKEFSA) $(LINK_EXTRA)   

GPGENMULT2 = \
  gpgenmult2.$O \
  $(LIBS)

$(BIN)/gpgenmult2: $(GPGENMULT2) 
	$(LINKER) -o $@ $(GPGENMULT2) $(LINK_EXTRA)   

GPMIGENMULT2 = \
  gpmigenmult2.$O \
  $(LIBS)

$(BIN)/gpmigenmult2: $(GPMIGENMULT2) 
	$(LINKER) -o $@ $(GPMIGENMULT2) $(LINK_EXTRA)   

GPGEOWA = \
  gpgeowa.$O \
  $(LIBS)

$(BIN)/gpgeowa: $(GPGEOWA) 
	$(LINKER) -o $@ $(GPGEOWA) $(LINK_EXTRA)   

GPORDER = \
  gporder.$O \
  $(LIBS)

$(BIN)/gporder: $(GPORDER) 
	$(LINKER) -o $@ $(GPORDER) $(LINK_EXTRA)   

GPOVLWA = \
  gpovlwa.$O \
  $(LIBS)

$(BIN)/gpovlwa: $(GPOVLWA) 
	$(LINKER) -o $@ $(GPOVLWA) $(LINK_EXTRA)   

GPSUBPRES = \
  gpsubpres.$O \
  $(LIBS)

$(BIN)/gpsubpres: $(GPSUBPRES) 
	$(LINKER) -o $@ $(GPSUBPRES) $(LINK_EXTRA)   

GPSTABILISER = \
  gpstabiliser.$O \
  $(LIBS)

$(BIN)/gpstabiliser: $(GPSTABILISER) 
	$(LINKER) -o $@ $(GPSTABILISER) $(LINK_EXTRA)


GPSUBLOWINDEX = \
  gpsublowindex.$O \
  $(LIBS)

$(BIN)/gpsublowindex: $(GPSUBLOWINDEX) 
	$(LINKER) -o $@ $(GPSUBLOWINDEX) $(LINK_EXTRA)

GPSUBMAKE = \
  gpsubmake.$O \
  $(LIBS)

$(BIN)/gpsubmake: $(GPSUBMAKE) 
	$(LINKER) -o $@ $(GPSUBMAKE) $(LINK_EXTRA)


GPSUBWA = \
  gpsubwa.$O \
  $(LIBS)

GPTCENUM = \
  gptcenum.$O \
  $(LIBS)

$(BIN)/gptcenum: $(GPTCENUM) 
	$(LINKER) -o $@ $(GPTCENUM) $(LINK_EXTRA) 

$(BIN)/gpsubwa: $(GPSUBWA) 
	$(LINKER) -o $@ $(GPSUBWA) $(LINK_EXTRA)   

GPVITAL = \
  gpvital.$O \
  $(LIBS)

$(BIN)/gpvital: $(GPVITAL) 
	$(LINKER) -o $@ $(GPVITAL) $(LINK_EXTRA)   

GPXLATWA = \
  gpxlatwa.$O \
  $(LIBS)

$(BIN)/gpxlatwa: $(GPXLATWA) 
	$(LINKER) -o $@ $(GPXLATWA) $(LINK_EXTRA)   

GPWA = \
  gpwa.$O \
  $(LIBS)

$(BIN)/gpwa: $(GPWA) 
	$(LINKER) -o $@ $(GPWA) $(LINK_EXTRA)   

GPDIFFLABS = \
  gpdifflabs.$O \
  $(LIBS)

$(BIN)/gpdifflabs: $(GPDIFFLABS) 
	$(LINKER) -o $@ $(GPDIFFLABS) $(LINK_EXTRA)   

GPMORPHISM = \
  gpmorphism.$O \
  $(LIBS)

$(BIN)/gpmorphism: $(GPMORPHISM) 
	$(LINKER) -o $@ $(GPMORPHISM) $(LINK_EXTRA)

FSALEQUAL = \
  fsalequal.$O \
  $(LIBS)

$(BIN)/fsalequal: $(FSALEQUAL) 
	$(LINKER) -o $@ $(FSALEQUAL) $(LINK_EXTRA)   

FSACARTESIAN = \
  fsacartesian.$O \
  $(LIBS)

$(BIN)/fsacartesian: $(FSACARTESIAN) 
	$(LINKER) -o $@ $(FSACARTESIAN) $(LINK_EXTRA)   

FSACOMPOSE = \
  fsacompose.$O \
  $(LIBS)

$(BIN)/fsacompose: $(FSACOMPOSE) 
	$(LINKER) -o $@ $(FSACOMPOSE) $(LINK_EXTRA)   

ISNORMAL = \
  isnormal.$O \
  $(LIBS)

$(BIN)/isnormal: $(ISNORMAL) 
	$(LINKER) -o $@ $(ISNORMAL) $(LINK_EXTRA)   

ISCONJUGATE = \
  isconjugate.$O \
  $(LIBS)

$(BIN)/isconjugate: $(ISCONJUGATE) 
	$(LINKER) -o $@ $(ISCONJUGATE) $(LINK_EXTRA)

SIMPLIFY = \
  simplify.$O \
  $(LIBS)

$(BIN)/simplify: $(SIMPLIFY) 
	$(LINKER) -o $@ $(SIMPLIFY) $(LINK_EXTRA)   

include mafu.dep
