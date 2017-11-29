#-------------------------------------------------
#
# Project created by QtCreator 2016-07-07T12:58:25
#
#-------------------------------------------------

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

greaterThan(QT_MAJOR_VERSION, 4){    
CONFIG += c++11
} else {
QMAKE_CXXFLAGS += -std=c++0x
}

TARGET = gp_builder

SOURCES += main.cpp\
    glycosylationsite.cpp \
    io.cpp \
    resolve_overlaps.cpp

HEADERS  += \
    glycosylationsite.h \
    io.h \
    resolve_overlaps.h

LIBS += -L$(GEMSHOME)/gmml/bin/ -lgmml

INCLUDEPATH += $(GEMSHOME)/gmml/includes/
INCLUDEPATH += $(GEMSHOME)/gmml/bin
DEPENDPATH += $(GEMSHOME)/gmml/bin
