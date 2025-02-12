
#ifndef REGIONPAIRDIALOG_H
#define REGIONPAIRDIALOG_H

#include "config/configuration.h"
#include "mesh/mesh.h"

#include <QDialog>

namespace espreso
{

namespace Ui {
class RegionPairDialog;
}

class RegionPairDialog : public QDialog
{
    Q_OBJECT

public:
    static RegionPairDialog* createRegionMaterial(ECFObject* map, Mesh* mesh, ECFObject* materials, ECFParameter* pair = nullptr);
    static RegionPairDialog* createRegionExpression(ECFObject* map, Mesh* mesh, ECFParameter* pair = nullptr);
    static RegionPairDialog* createRegionOption(ECFObject* map, Mesh* mesh, ECFParameter* pair = nullptr);
    static RegionPairDialog* createRegionObject(ECFObject* map, Mesh* mesh, ECFObject* pair = nullptr);
    ~RegionPairDialog();

    void accept() override;
    void reject() override;

    QString region();

private:
    explicit RegionPairDialog(ECFDataType value,
                             ECFObject* map, Mesh* mesh,
                             ECFObject* scope, QWidget *parent = 0);
    explicit RegionPairDialog(ECFParameter* pair, ECFDataType value,
                             ECFObject* map, Mesh* mesh,
                             ECFObject* scope, QWidget *parent = 0);

    explicit RegionPairDialog(ECFObject* map, Mesh* mesh,
                              QWidget* parent = 0);
    explicit RegionPairDialog(ECFObject* pair, ECFObject* map,
                              Mesh* mesh, QWidget* parent = 0);

    Ui::RegionPairDialog *ui;

    Mesh* m_mesh;
    ECFObject* m_scope;
    ECFObject* m_map;
    ECFDataType m_first;
    ECFDataType m_second;
    QWidget* m_first_widget;
    QWidget* m_second_widget;

    ECFObject* m_object = nullptr;

    QString m_region;

    QWidget* uiValue(ECFDataType type);

    void displayError(const QString&);

    std::string getKey();
};

}

#endif // REGIONPAIRDIALOG_H
