#ifndef DATATYPEEDITWIDGET_H
#define DATATYPEEDITWIDGET_H

#include <QWidget>

#include <QLineEdit>
#include <QComboBox>
#include "tabletypewidget.h"
#include "piecewisetypewidget.h"
#include "config/configuration.h"
#include "elements/expressionedit.h"
#include "elements/ivalidatableobject.h"
#include "elements/isavableobject.h"
#include "elements/textitemwidget.h"

namespace espreso
{

    namespace Ui {
    class DataTypeEditWidget;
    }

	struct DataTypeEditWidgetFactoryData
	{
		int type;
		bool valid;
		QString error_message;
	};

    class DataTypeEditWidget : public TextItemWidget, public IValidatableObject,
            public ISavableObject
    {
        Q_OBJECT

    public:
        static QStringList typeNames();

        explicit DataTypeEditWidget(ECFParameter* data, QWidget* parent = 0);
        explicit DataTypeEditWidget(const std::vector<std::string>& variables, QWidget* parent = 0);
        ~DataTypeEditWidget();

        QComboBox* createComboBox(QWidget* parent = nullptr);
        void setComboBox(bool show);

        bool isValid() override;
        QString errorMessage() override;
        void save() override;

        virtual void setText(const QString& text) override;
        virtual QString text() override;

        QString value();

		void setSharedData(DataTypeEditWidgetFactoryData* data);

		int datatype();

    private slots:
        void changeType(int index);
		void onCellChange();

    private:
        DataTypeEditWidget(QWidget *parent = 0);

        Ui::DataTypeEditWidget *ui;

        ECFParameter* m_param;
		QString m_text;

        QComboBox* m_cmb;
        QComboBox* m_external_cmb = nullptr;

        ExpressionEdit* uiExpression;
        TableTypeWidget* uiTable;
        PiecewiseTypeWidget* uiPiecewise;

		DataTypeEditWidgetFactoryData* m_shared = nullptr;

		QVector<bool> m_datatype_empty;

        int activeType;
        std::vector<std::string> m_param_variables;
        QString param_getValue();
        std::vector<std::string> param_variables();

		void createUi();
		void parseValue(const QString&);
		int checkDataType(const QString&);
        void initExpression();
        void initTable();
        void initPiecewise();
    };

}

#endif // DATATYPEEDITWIDGET_H
