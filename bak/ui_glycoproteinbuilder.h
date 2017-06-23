/********************************************************************************
** Form generated from reading UI file 'glycoproteinbuilder.ui'
**
** Created by: Qt User Interface Compiler version 5.8.0
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_GLYCOPROTEINBUILDER_H
#define UI_GLYCOPROTEINBUILDER_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_glycoproteinBuilder
{
public:
    QMenuBar *menuBar;
    QToolBar *mainToolBar;
    QWidget *centralWidget;
    QStatusBar *statusBar;

    void setupUi(QMainWindow *glycoproteinBuilder)
    {
        if (glycoproteinBuilder->objectName().isEmpty())
            glycoproteinBuilder->setObjectName(QStringLiteral("glycoproteinBuilder"));
        glycoproteinBuilder->resize(400, 300);
        menuBar = new QMenuBar(glycoproteinBuilder);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        glycoproteinBuilder->setMenuBar(menuBar);
        mainToolBar = new QToolBar(glycoproteinBuilder);
        mainToolBar->setObjectName(QStringLiteral("mainToolBar"));
        glycoproteinBuilder->addToolBar(mainToolBar);
        centralWidget = new QWidget(glycoproteinBuilder);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        glycoproteinBuilder->setCentralWidget(centralWidget);
        statusBar = new QStatusBar(glycoproteinBuilder);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        glycoproteinBuilder->setStatusBar(statusBar);

        retranslateUi(glycoproteinBuilder);

        QMetaObject::connectSlotsByName(glycoproteinBuilder);
    } // setupUi

    void retranslateUi(QMainWindow *glycoproteinBuilder)
    {
        glycoproteinBuilder->setWindowTitle(QApplication::translate("glycoproteinBuilder", "glycoproteinBuilder", Q_NULLPTR));
    } // retranslateUi

};

namespace Ui {
    class glycoproteinBuilder: public Ui_glycoproteinBuilder {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_GLYCOPROTEINBUILDER_H
