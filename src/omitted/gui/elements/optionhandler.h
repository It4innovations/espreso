#ifndef OPTIONHANDLER_H
#define OPTIONHANDLER_H

#include <QWidget>

#include "config/configuration.h"

namespace espreso {
    namespace Ui {
        class OptionHandler;
    }

    class OptionHandler : public QWidget
    {
        Q_OBJECT

    public:
        explicit OptionHandler(ECFParameter*, QWidget*, bool = true);
        ~OptionHandler();

    signals:
        void optionChanged();

    private slots:
        void onIndexChanged(int index);

    private:
        Ui::OptionHandler *ui;
        ECFParameter* m_option;
        bool optionsAdded = false;
    };
}

#endif // OPTIONHANDLER_H
