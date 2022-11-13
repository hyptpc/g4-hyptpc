#
#	Makefile
#
# Compilers
C++	= g++
CC	= gcc
#
# geant4.9.3 #geant4.9.6
GEANT4DIR	= /home/sks/software/geant4/geant4.9.6.p02
#
#G4LIBDIR	= $(GEANT4DIR)/lib/Linux-g++
#G4INCDIR	= $(GEANT4DIR)/include
G4LIBDIR	= $(GEANT4DIR)/lib64/Geant4-9.6.2/Linux-g++
G4INCDIR	= $(GEANT4DIR)/include/Geant4

#
G4INCLUDE	= -I$(G4INCDIR)
G4LIBS		= -L$(G4LIBDIR)			\
			-lG4global		\
			-lG4digits_hits		\
			-lG4geometry		\
			-lG4interfaces		\
			-lG4parmodels		\
			-lG4particles		\
			-lG4processes		\
			-lG4event		\
			-lG4graphics_reps	\
			-lG4intercoms		\
			-lG4materials		\
			-lG4modeling		\
			-lG4readout		\
			-lG4run			\
			-lG4track		\
			-lG4tracking		\
			-lG4vis_management	\
			-lG4FR			\
			-lG4VRML		\
			-lG4OpenGL		\
			-lG4RayTracer		\
			-lG4visXXX		\
			-lG4gl2ps		
#                       -lG4emstandard          
#
CLHEPDIR	= /home/sks/software/CLHEP/clhep-2.1.3.1
CLHEPINCDIR	= $(CLHEPDIR)/include/CLHEP
CLHEPLIBDIR	= $(CLHEPDIR)/lib
CLHEPINCLUDE	= -I$(CLHEPINCDIR)
CLHEPLIB	= -L$(CLHEPLIBDIR) -lCLHEP
#
#
EVTGENDIR	= /home/sks/software/EvtGen/EvtGen/EvtGen/R01-01-00
EVTGENINCDIR	= $(EVTGENDIR)
EVTGENLIBDIR	= $(EVTGENDIR)/lib
EVTGENINCLUDE	= -I$(EVTGENINCDIR)
EVTGENLIB	= -L$(EVTGENLIBDIR) -lEvtGen -lEvtGenExternal
EVTGEN_EXTERNAL = 0
EXTRAFLAGS      = -DEVTGEN_EXTERNAL
#
#
HEPMCDIR	= /usr/local
HEPMCINCDIR	= $(HEPMCDIR)/include/HepMC
HEPMCLIBDIR	= $(HEPMCDIR)/lib
HEPMCINCLUDE	= -I$(HEPMCINCDIR)
HEPMCLIB	= -L$(HEPMCLIBDIR) -lHepMC
#
#
PYTHIADIR	= /home/sks/software/EvtGen/PYTHIA/pythia8170
PYTHIAINCDIR	= $(PYTHIADIR)/include
PYTHIALIBDIR	= $(PYTHIADIR)/lib
PYTHIAINCLUDE	= -I$(PYTHIAINCDIR)
PYTHIALIB	= -L$(PYTHIALIBDIR) -lpythia8 -llhapdfdummy -lhepmcinterface
#
#
PHOTOSDIR	= /home/sks/software/EvtGen/PHOTOS/PHOTOS
PHOTOSINCDIR	= $(PHOTOSDIR)/include
PHOTOSLIBDIR	= $(PHOTOSDIR)/lib
PHOTOSINCLUDE	= -I$(PHOTOSINCDIR)
PHOTOSLIB	= -L$(PHOTOSLIBDIR) -lPhotosCxxInterface -lPhotosFortran
#
#
#X11DIR		= /usr/X11R6
X11DIR		= /usr
X11LIBDIR	= $(X11DIR)/lib
X11INCDIR	= $(X11DIR)/include
#
XINCLUDE	= -I$(X11INCDIR)
XLIBS		= -L$(X11LIBDIR) -lX11 -lXmu
#
#MESADIR		= /usr/X11R6/lib
MESADIR		= /usr/lib
MESALIBS	= -L$(MESADIR) -lGL -lGLU
#
#
BOOSTDIR	= /usr
BOOSTLIBDIR	= $(BOOSTDIR)/lib64
BOOSTINCDIR	= $(BOOSTDIR)/include
BOOSTINCLUDE	= -I$(BOOSTINCDIR)
BOOSTLIB	= -L$(BOOSTLIBDIR) -lboost_iostreams
#
#
RTINCLUDE	= `root-config --cflags`
RTLIBS		= `root-config --libs`
#
DEFINEFLAGS	=	-DG4UI_USE_TCSH		\
			-DG4UI_USE_TERMINAL	\
			-DG4VIS_USE		\
			-DG4VIS_USE_OPENGL	\
			-DG4VIS_USE_OPENGLX	\
			-DG4VIS_USE_VRMLFILE	\
			-DG4VIS_USE_DAWNFILE	\
			-DG4VIS_USE_RAYTRACER	\
			-DG4USE_STD_NAMESPACE	\
			-DEVTGEN_EXTERNAL
#
DEBUGFLAGS	= -g -Wall
#
CPPFLAGS	= $(DEFINEFLAGS) $(DEBUGFLAGS) -I. $(G4INCLUDE)		\
			$(CLHEPINCLUDE) $(XINCLUDE) $(RTINCLUDE)	\
			$(EVTGENINCLUDE) $(HEPMCINCLUDE)		\
			$(PYTHIAINCLUDE) $(PHOTOSINCLUDE)		\
			$(BOOSTINCLUDE)
#
LIBS		= $(G4LIBS) $(CLHEPLIB) $(MESALIBS) $(XLIBS)		\
			$(CERNLIBS) $(RTLIBS)	$(EVTGENLIB) $(HEPMCLIB)	\
			$(PYTHIALIB) $(PHOTOSLIB) $(BOOSTLIB)
#
CFLAGS		= -O2
#
FLAGS		= $(CFLAGS) $(CPPFLAGS)
#
BINDIR		= ../bin
#
TOUCH		= touch
#
RM		= rm -f
#
# Rules
#
.cc.o:
	$(C++) $(FLAGS) -MMD -c $<
.cpp.o:
	$(C++) $(FLAGS) -MMD -c $<
#
.c.o:
	$(CC) $(FLAGS) -MMD -c $<
#
#

OBJS	= 	GetNumberFromKernelEntropyPool.o	\
		ConfMan.o				\
		SpecField.o				\
		SpecFieldMap.o				\
		SimpleFieldElement.o			\
		DetectorConstruction.o			\
		PhysicsList.o				\
		RunAction.o				\
		PrimaryGeneratorAction.o		\
		GeneratorHelper.o			\
		SteppingAction.o			\
		Analysis.o				\
		Particle.o				\
		VisManager.o				\
		Transportation.o			

all: $(BINDIR)/main 

clean:
	$(RM) $(OBJS) $(OBJS_SPEC) *.d *~
#
SRCS	= $(wildcard *.cc)
-include $(SRCS:%.cc=%.d)

$(BINDIR)/main: $(OBJS)
	$(C++) -o $@ $(CFLAGS) $(OBJS) $(LIBS)

