#-------------------------------------------------
#
# Project created by QtCreator 2016-07-07T12:58:25
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

greaterThan(QT_MAJOR_VERSION, 4){    
CONFIG += c++11
} else {
QMAKE_CXXFLAGS += -std=c++0x
}

TARGET = glycoproteinBuilder
TEMPLATE = app


SOURCES += main.cpp\
    glycosylationsite.cpp \
    overlap.cpp \
    io.cpp \
    attachedrotamer.cpp

HEADERS  += \
    glycosylationsite.h \
    overlap.h \
    io.h \
    attachedrotamer.h

FORMS    += glycoproteinbuilder.ui




LIBS += -L/home/oliver/Programs/gems/gmml/bin/ -lgmml

INCLUDEPATH += $$PWD/../../../gems/gmml/includes/
INCLUDEPATH += $$PWD/../../../gems/gmml/bin
DEPENDPATH += $$PWD/../../../gems/gmml/bin
