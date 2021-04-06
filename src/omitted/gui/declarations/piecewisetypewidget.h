#ifndef PIECEWISETYPEWIDGET_H
#define PIECEWISETYPEWIDGET_H

#include "tablewidget.h"

namespace espreso
{

    class PiecewiseTypeWidget : public TableWidget
    {
        Q_OBJECT

    public:
        static QStringList headlines();

        PiecewiseTypeWidget(const std::vector<std::string>& variables,
                            QWidget* parent = 0);
        bool isValid() override;
        QString errorMessage() override;
        QString data() override;

        virtual void addData(const QString&) override;

    private:
        std::vector<std::string> m_variables;
        int m_invalidRow = 0;
        bool m_isValid = true;
    };

}

#endif // PIECEWISETYPEWIDGET_H
