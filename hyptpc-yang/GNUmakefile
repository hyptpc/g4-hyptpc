name := main
G4TARGET := $(name)
G4EXLIB := true

.PHONY: all
all:  lib bin

include $(G4INSTALL)/config/binmake.gmk

#root
CPPFLAGS += `root-config --cflags`
LDFLAGS  += `root-config --libs`

#ROOTLIBS       = $(shell root-config --libs)
ROOTDIR     =   /software/ROOT/ROOT-v6.24.06/root-6.24.06-install
ROOTINCDIR  =   $(ROOTDIR)/include
ROOTLIBDIR  =   $(ROOTDIR)/lib
ROOTINCLUDE =  -I$(ROOTINCDIR)
ROOTLIB =  -L$(ROOTLIBDIR)


CPPFLAGS += $(ROOTINCLUDE)
LDFLAGS  += $(ROOTLIB)
#LDFLAGS  += $(ROOTLIBS)

CPPFLAGS += $(shell $(ROOTDIR)/bin/root-config --cflags)
EXTRALIBS += $(shell $(ROOTDIR)/bin/root-config --libs)
LDLIBS += $(shell $(ROOTDIR)/bin/root-config --libs)


#event gen
EVTGENDIR       = /home/E72/software/EvtGen/evtgen-install
EVTGENINCDIR    = $(EVTGENDIR)/include
EVTGENLIBDIR    = $(EVTGENDIR)/lib
EVTGENINCLUDE   = -I$(EVTGENINCDIR)
EVTGENLIB       = -L$(EVTGENLIBDIR) -lEvtGen -lEvtGenExternal
EVTGEN_EXTERNAL = 0
EXTRAFLAGS      = -DEVTGEN_EXTERNAL

CPPFLAGS += $(EVTGENINCLUDE)
LDFLAGS += $(EVTGENLIB)
LDLIBS += $(EVTGENLIB)

#HEPMC
HEPMCDIR       = /home/E72/software/HepMC/HepMC-2.06.09-install
HEPMCINCDIR    = $(HEPMCDIR)/include
HEPMCLIBDIR    = $(HEPMCDIR)/lib
HEPMCINCLUDE   = -I$(HEPMCINCDIR)
HEPMCLIB       = -L$(HEPMCLIBDIR) -lHepMC

CPPFLAGS += $(HEPMCINCLUDE)
LDFLAGS += $(HEPMCLIB)
LDLIBS += $(HEPMCLIB)

#Pythia
PYTHIADIR       = /home/E72/software/PYTHIA/pythia8245
PYTHIAINCDIR    = $(PYTHIADIR)/include
PYTHIALIBDIR    = $(PYTHIADIR)/lib
PYTHIAINCLUDE   = -I$(PYTHIAINCDIR)
#PYTHIALIB       = -L$(PYTHIALIBDIR) -lpythia8 -llhapdfdummy -lhepmcinterface
PYTHIALIB       = -L$(PYTHIALIBDIR) -lpythia8

CPPFLAGS += $(PYTHIAINCLUDE)
LDFLAGS += $(PYTHIALIB)
LDLIBS += $(PYTHIALIB)

#PHOTOS                                                                                                     
PHOTOSDIR       = /home/E72/software/PHOTOS/PHOTOS
PHOTOSINCDIR    = $(PHOTOSDIR)/include
PHOTOSLIBDIR    = $(PHOTOSDIR)/lib
PHOTOSINCLUDE   = -I$(PHOTOSINCDIR)
PHOTOSLIB       = -L$(PHOTOSLIBDIR) -lPhotosCxxInterface -lPhotosFortran

CPPFLAGS += $(PHOTOSINCLUDE)
LDFLAGS += $(PHOTOSLIB)
LDLIBS += $(PHOTOSLIB)


#CLHEP
CLHEPDIR        = /software/CLHEP/CLHEP2.4.5.1/CLHEP_install
CLHEPINCDIR     = $(CLHEPDIR)/include/CLHEP
CLHEPLIBDIR     = $(CLHEPDIR)/lib
CLHEPINCLUDE    = -I$(CLHEPINCDIR)
CLHEPLIB        = -L$(CLHEPLIBDIR) -lCLHEP

CPPFLAGS += $(CLHEPINCLUDE)
LDFLAGS += $(CLHEPLIB)
LDLIBS += $(CLHEPLIB)

#X11
#X11DIR         = /usr/X11R6                                                                      
X11DIR          = /usr
X11LIBDIR       = $(X11DIR)/lib64
X11INCDIR       = $(X11DIR)/include
XINCLUDE        = -I$(X11INCDIR)
XLIBS           = -L$(X11LIBDIR) -lX11 -lXmu

CPPFLAGS += $(XINCLUDE)
LDFLAGS += $(XLIB)

CPPFLAGS += -DG4MULTITHREADED
