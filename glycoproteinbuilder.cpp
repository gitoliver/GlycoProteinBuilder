
#include "glycoproteinbuilder.h"
#include "ui_glycoproteinbuilder.h"

glycoproteinBuilder::glycoproteinBuilder(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::glycoproteinBuilder)
{
    ui->setupUi(this);
}

glycoproteinBuilder::~glycoproteinBuilder()
{
    delete ui;
}
