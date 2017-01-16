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
        glycoproteinbuilder.cpp \
    gpb_inputs.cpp

HEADERS  += glycoproteinbuilder.h \
    gpb_inputs.h

FORMS    += glycoproteinbuilder.ui

win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../gems/gmml/bin/release/ -lgmml
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../gems/gmml/bin/debug/ -lgmml
else:unix: LIBS += -L$$PWD/../../../gems/gmml/bin/ -lgmml


LIBS += -L/home/oliver/Programs/gems/gmml/bin/ -lgmml

INCLUDEPATH += /home/oliver/Programs/gems-Nov2016/gmml/includes/Eigen_Algebra_Templates
INCLUDEPATH += $$PWD/../../../gems/gmml/bin
DEPENDPATH += $$PWD/../../../gems/gmml/bin
