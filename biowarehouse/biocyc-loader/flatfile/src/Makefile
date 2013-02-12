all: biocyc-loader
db := none
utilsdir := ../../../utils/src/c
VPATH = $(utilsdir)

## If we don't specify a $(db), make sure both mysql and oracle files are clean'ed
OBJS = wh_oracle.o wh_oracle_util.o wh_mysql.o wh_mysql_util.o

ifeq ($(db),oracle)
	##
	## Include Oracle make rules
	##
	include $(ORACLE_HOME)/precomp/lib/env_precomp.mk
	MODE=oracle
	HDRS = wh_oracle.h wh_oracle_util.h
        OBJS = wh_oracle.o wh_oracle_util.o
        LIBS = $(LDPATHFLAG)$(LIBHOME) $(PROLDLIBS)
        OUTPUT = oracle-biocyc-loader
        CFLAGS += -I$(utilsdir)
endif
ifeq ($(db), mysql)
        HDRS = wh_mysql.h
        OBJS = wh_mysql.o wh_mysql_util.o
        LIBS = -L/usr/lib -lmysqlclient -lz -lm -lcrypt -lpthread -ldl
        OUTPUT = mysql-biocyc-loader
        CFLAGS += -I/usr/include/mysql -I$(utilsdir)
endif

CC=gcc
HDRS +=main.h db.h string-util.h widtable.h wh.h 
OBJS +=main.o string-util.o widtable.o db.o wh.o\
	pubs.tab.o pubs-load.o \
	compound.tab.o compound-load.o \
	reaction.tab.o reaction-load.o \
	protein.tab.o protein-load.o \
	protseq.tab.o protseq-load.o \
	enzrxn.tab.o enzrxn-load.o \
	regulation.tab.o regulation-load.o \
	gene.tab.o gene-load.o \
	transunit.tab.o transunit-load.o \
	promoter.tab.o promoter-load.o \
	terminator.tab.o terminator-load.o \
	dnabindsite.tab.o dnabindsite-load.o \
	pathway.tab.o pathway-load.o 

## For gdb 
CFLAGS+= -g

biocyc-loader: $(OBJS) $(HDRS)
	$(CC) -O3 -g -o $(OUTPUT) $(OBJS) -ll $(LIBS) ## -lefence
	rm -f biocyc-loader
	ln -s $(OUTPUT) biocyc-loader

##
## Pubs
##
pubs.tab.c:	pubs.y lex.pubs.c
	bison -v -p pubs pubs.y

lex.pubs.c: pubs.l
	flex -Ppubs -s pubs.l

pubs.tab.o: lex.pubs.c $(HDRS)
pubs-load.o: pubs-parse.h

##
## Compound
##
compound.tab.c:	compound.y lex.compound.c
	bison -v -p compound compound.y

lex.compound.c: compound.l
	flex -Pcompound -s compound.l

compound.tab.o: lex.compound.c $(HDRS)
compound-load.o: compound-parse.h

##
## Reaction
##
reaction.tab.c:	reaction.y lex.reaction.c
	bison -v -p reaction reaction.y

lex.reaction.c: reaction.l
	flex -Preaction -s reaction.l

reaction.tab.o: lex.reaction.c $(HDRS)
reaction-load.o: reaction-parse.h

##
## Pathway
##
pathway.tab.c:	pathway.y lex.pathway.c
	bison -v -p pathway pathway.y

lex.pathway.c: pathway.l
	flex -Ppathway -s pathway.l

pathway.tab.o: lex.pathway.c $(HDRS)
pathway-load.o: pathway-parse.h

##
## Protein
##
protein.tab.c:	protein.y lex.protein.c
	bison -v -p protein protein.y

lex.protein.c: protein.l
	flex -Pprotein -s protein.l

protein.tab.o: lex.protein.c $(HDRS)
protein-load.o: protein-parse.h

##
## Protseq
##
protseq.tab.c:	protseq.y lex.protseq.c
	bison -v -p protseq protseq.y

lex.protseq.c: protseq.l
	flex -Pprotseq -s protseq.l

protseq.tab.o: lex.protseq.c $(HDRS)
protseq-load.o: protseq-parse.h

## Gene
##
gene.tab.c:	gene.y lex.gene.c
	bison -v -p gene gene.y

lex.gene.c: gene.l
	flex -Pgene -s gene.l

gene.tab.o: lex.gene.c $(HDRS)
gene-load.o: gene-parse.h

## Transunit
##
transunit.tab.c:	transunit.y lex.transunit.c
	bison -v -p transunit transunit.y

lex.transunit.c: transunit.l
	flex -Ptransunit -s transunit.l

transunit.tab.o: lex.transunit.c $(HDRS)
transunit-load.o: transunit-parse.h

## Promoter
##
promoter.tab.c:	promoter.y lex.promoter.c
	bison -v -p promoter promoter.y

lex.promoter.c: promoter.l
	flex -Ppromoter -s promoter.l

promoter.tab.o: lex.promoter.c $(HDRS)
promoter-load.o: promoter-parse.h

## Terminator
##
terminator.tab.c:	terminator.y lex.terminator.c
	bison -v -p terminator terminator.y

lex.terminator.c: terminator.l
	flex -Pterminator -s terminator.l

terminator.tab.o: lex.terminator.c $(HDRS)
terminator-load.o: terminator-parse.h

## DNAbindsite
##
dnabindsite.tab.c:	dnabindsite.y lex.dnabindsite.c
	bison -v -p dnabindsite dnabindsite.y

lex.dnabindsite.c: dnabindsite.l
	flex -Pdnabindsite -s dnabindsite.l

dnabindsite.tab.o: lex.dnabindsite.c $(HDRS)
dnabindsite-load.o: dnabindsite-parse.h

##
## Enzrxn
##
enzrxn.tab.c: enzrxn.y lex.enzrxn.c
	bison -v -p enzrxn enzrxn.y

lex.enzrxn.c: enzrxn.l
	flex -Penzrxn -s enzrxn.l

enzrxn.tab.o: lex.enzrxn.c $(HDRS)
enzrxn-load.o: enzrxn-parse.h

## Regulation
##
regulation.tab.c: regulation.y lex.regulation.c
	bison -v -p regulation regulation.y

lex.regulation.c: regulation.l
	flex -Pregulation -s regulation.l

regulation.tab.o: lex.regulation.c $(HDRS)
regulation-load.o: regulation-parse.h

##
## Supporting Rules
##
.PHONY: clean 

clean:
	rm -f $(OBJS) $(utilsdir)/*.lis $(utilsdir)/wh_oracle_util.c wh_oracle.o wh_oracle_util.o wh_mysql.o wh_mysql_util.o compound-load.c pathway-load.c protein-load.c protseq-load.c enzrxn-load.c gene-load.c reaction-load.c promoter-load.c transunit-load.c dnabindsite-load.c biocyc-loader lex.*.c *.tab.c *.lis *.output $(OUTPUT)

##
## Compilation rules for esql programs
##
ifeq ($(db), oracle)
.SUFFIXES: .pc
.pc.c:
	$(PROC) DEFINE=ORACLE INCLUDE=$(utilsdir) $(PROCFLAGS) iname=$<
.c.o:
	$(CC) -DORACLE $(CFLAGS) -c $<
endif
ifeq ($(db), mysql)
%.c:%.pc
	cp $*.pc $*.c
%.o:%.c
	$(CC) -DDEF_MYSQL $(CFLAGS) -c $<
endif
