#ifndef MATERIALDIALOG_H
#define MATERIALDIALOG_H

#include <QDialog>
#include <QFrame>
#include <QVBoxLayout>
#include <QVector>
#include "config/ecf/material/material.h"
#include "elements/isavableobject.h"
#include "elements/ivalidatableobject.h"
#include "materialwidget.h"

namespace espreso
{

    namespace Ui {
    class MaterialDialog;
    }

    class MaterialDialog : public QDialog
    {
        Q_OBJECT

    public:
        explicit MaterialDialog(MaterialConfiguration* material,
                                const QVector<std::string>& materialNames,
                                QWidget *parent = 0);
        ~MaterialDialog();

        void accept() override;

    private:
        Ui::MaterialDialog *ui;

        MaterialWidget* m_widget;
    };

}

#endif // MATERIALDIALOG_H
