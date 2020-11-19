#ifndef MATERIALDIALOG_H
#define MATERIALDIALOG_H

#include "materialwidget.h"
#include "gui/elements/isavableobject.h"
#include "gui/elements/ivalidatableobject.h"

#include "config/ecf/material/material.h"

#include <QDialog>
#include <QFrame>
#include <QVBoxLayout>
#include <QVector>

namespace espreso
{

    namespace Ui {
    class MaterialDialog;
    }

    class MaterialDialog : public QDialog
    {
        Q_OBJECT

    public:
        explicit MaterialDialog(ECFObject* material,
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
