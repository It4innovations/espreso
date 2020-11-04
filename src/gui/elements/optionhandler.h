
#ifndef OPTIONHANDLER_H
#define OPTIONHANDLER_H

#include "config/configuration.h"

#include <QWidget>

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
