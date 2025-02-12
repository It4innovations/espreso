
#include "textitemwidgetfactory.h"

using namespace espreso;

TextItemWidgetFactory::TextItemWidgetFactory()
{

}


TextItemWidget* FilepathWidgetFactory::create(QWidget *parent)
{
    return new FilepathWidget(parent);
}


OptionWidgetFactory::OptionWidgetFactory(const QStringList &options)
{
    this->m_options = options;
}

TextItemWidget* OptionWidgetFactory::create(QWidget *parent)
{
    return new OptionWidget(this->m_options, parent);
}

TextWidgetFactory::TextWidgetFactory(ValidatorFactory *validatorFactory)
{
    this->m_factory = validatorFactory;
}

TextWidgetFactory::~TextWidgetFactory()
{
    if (this->m_factory != nullptr) delete this->m_factory;
}

TextItemWidget* TextWidgetFactory::create(QWidget *parent)
{
    TextWidget* widget = new TextWidget(parent);
    if (this->m_factory != nullptr)
        widget->setValidator(this->m_factory->create());
    return widget;
}

DataTypeEditWidgetFactory::DataTypeEditWidgetFactory(ECFParameter* expression)
{
    this->m_expr = expression;
    if (m_expr->getValue().empty()) this->m_expr->setValue("0");

    this->m_data.valid = true;

    DataTypeEditWidget* w = new DataTypeEditWidget(this->m_expr);
    int datatype = w->datatype();
    w->hide();
    delete w;
    this->m_data.type = datatype;
}

TextItemWidget* DataTypeEditWidgetFactory::create(QWidget* parent)
{
    DataTypeEditWidget* widget = new DataTypeEditWidget(this->m_expr->metadata.variables, parent);
    widget->setComboBox(true);
    widget->setSharedData(&this->m_data);

    return widget;
}

bool DataTypeEditWidgetFactory::isValid()
{
    return this->m_data.valid;
}

QString DataTypeEditWidgetFactory::errorMessage()
{
    return this->m_data.error_message;
}
