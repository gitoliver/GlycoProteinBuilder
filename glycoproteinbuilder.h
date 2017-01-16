
#ifndef GLYCOPROTEINBUILDER_H
#define GLYCOPROTEINBUILDER_H

#include <QMainWindow>
namespace Ui {
class glycoproteinBuilder;
}

class glycoproteinBuilder : public QMainWindow
{
    Q_OBJECT

public:
    explicit glycoproteinBuilder(QWidget *parent = 0);
    ~glycoproteinBuilder();

private:
    Ui::glycoproteinBuilder *ui;
};

#endif // GLYCOPROTEINBUILDER_H
